#!/usr/bin/env python3
'''
Preprocessing script running the preprocessing based on configuration
options.
'''
# standard libraries
import os
import logging
import logging.config
from typing import Tuple, List

import pickle
from argparse import ArgumentParser
import json

import yaml

# own libraries
from lib import download, errorfixer
from lib.visual import progress_bar
from lib.model import json_parser, case, config
from lib.api import phenomizer, omim, mutalyzer


def configure_logging(logger_name, logger_file: str = "preprocess.log"):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    # visible screen printing
    stdout_channel = logging.StreamHandler()
    stdout_channel.setLevel(logging.INFO)
    # file output logging
    file_channel = logging.FileHandler(logger_file, mode="w")
    file_channel.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(message)s")
    stdout_channel.setFormatter(formatter)

    file_formatter = logging.Formatter(
        "%(asctime)s L%(lineno)d <%(module)s|%(funcName)s> %(message)s"
    )
    file_channel.setFormatter(file_formatter)

    logger.addHandler(stdout_channel)
    logger.addHandler(file_channel)


def parse_arguments():
    parser = ArgumentParser(description=(
        "Process f2g provided jsons into a format processable by "
        "classification."))
    parser.add_argument("-s", "--single", help="Process a single json file.")
    parser.add_argument(
        "-o", "--output",
        help="Destination of created old json.",
        default=""
    )
    parser.add_argument(
        "-p", "--pickle",
        help="Start with pickled cases after phenomization."
    )
    parser.add_argument(
        "-e", "--entry",
        help=("Start entrypoint for pickled results. "
              "Default: pheno - start at phenomization. "
              "convert - start at old json mapping. "
              "qc - start at case quality check. "
              "Used in conjunction with --pickle."),
        default="pheno"
    )
    parser.add_argument(
        "--skip-vcf", action='store_true',
        help="Skip vcf convertion."
    )

    return parser.parse_args()


def json_from_directory(config_data: config.ConfigManager) \
        -> Tuple[List[str], str]:
    '''Get a list of json file paths.'''
    # Download new files from AWS Bucket
    if config_data.general.getboolean("download"):
        files = download.backup_s3_folder(config=config_data)
    else:
        files = None

    # Initial Quality check of new json
    unprocessed_jsons = os.path.join(
        config_data.aws['download_location'], 'cases'
    )
    if config_data.general.getboolean("changed_only") and files is not None:
        json_files = [
            os.path.join(unprocessed_jsons, f)
            for f in files["cases"]
        ]
    else:
        json_files = [
            os.path.join(unprocessed_jsons, x)
            for x in os.listdir(unprocessed_jsons)
            if os.path.splitext(x)[1] == '.json'
        ]
    # corrected is a directory which can contain manually edited case jsons
    # that should differ from the original only in content, not in overall
    # structure.
    # this should make resolving some exotic errors a lot easier
    corrected = config_data.preprocess['corrected_location']

    return json_files, corrected


def create_config(
        config_path: str = "config.yml",
        simvcffolder: str = "data/PEDIA/mutations",
        vcffolder: str = "data/PEDIA/vcfs/original",
        merge: bool = True
) -> None:
    '''Creates config.yml file based on the VCF files'''
    # real vcf files
    vcffiles = {
        f.split(".")[0]
        for f in os.listdir(vcffolder) if not f.startswith(".")
    }
    # simulated vcf files
    singlefiles = {
        f.split(".")[0]
        for f in os.listdir(simvcffolder) if not f.startswith(".")
    }

    # completely simulated cases
    single_samples = list(singlefiles - vcffiles)
    # cases with vcfs and simulated data
    vcf_samples = list(vcffiles & singlefiles)
    # cases only with vcfs are used for testing
    test_samples = list(vcffiles - singlefiles)

    config_data = {
        "SINGLE_SAMPLES": single_samples,
        "VCF_SAMPLES": vcf_samples,
        "TEST_SAMPLES": test_samples
    }
    if merge and os.path.exists(config_path):
        with open(config_path, "r") as configfile:
            olddata = yaml.load(configfile)
        config_data = {
            k: list(set(v) | set(olddata[k]))
            for k, v in config_data.items()
        }

    with open("config.yml", "w") as configfile:
        yaml.dump(config_data, configfile, default_flow_style=False)


def create_jsons(args, config_data):
    '''Create a list of new formatjson objects.'''
    print("== Process new json files ==")
    # get either from single file or from directory
    json_files, corrected = ([args.single], config_data.preprocess['corrected_location']) \
        if args.single else json_from_directory(config_data)

    @progress_bar("Process jsons")
    def yield_jsons(json_files, corrected):
        for json_file in json_files:
            yield json_parser.NewJson.from_file(json_file, corrected)

    new_json_objs = yield_jsons(json_files, corrected)

    print('Unfiltered', len(new_json_objs))

    if config_data.jsonparser["json_qc_log"] and not args.single:
        logpath = config_data.jsonparser["json_qc_log"]
        qc_failed_results = {
            case_id: {
                "issues": issues,
                "valid": valid
            }
            for case_id, (valid, issues) in
            [(j.get_case_id(), j.check()) for j in new_json_objs]
        }
        # merge with old log data if partial update
        if config_data.general.getboolean("changed_only") \
                and os.path.exists(logpath):
            with open(logpath, "r") as failedfile:
                olddata = json.load(failedfile)
            qc_failed_results = {**olddata, **qc_failed_results}
        with open(logpath, "w") as failedfile:
            json.dump(qc_failed_results, failedfile, indent=4)

    filtered_new = [j for j in new_json_objs if j.check()[0]]
    print('Filtered rough criteria', len(filtered_new))
    return filtered_new


def create_cases(args, config_data, jsons):
    '''Create cases from list of jsons.'''
    print("== Create cases from new json format ==")
    error_fixer = errorfixer.ErrorFixer(config=config_data)
    omim_obj = omim.Omim(config=config_data)

    @progress_bar("Create cases")
    def yield_cases(json_files, error_fixer, omim_obj, exclusion):
        '''Create case from json objects.'''
        for json_file in json_files:
            yield case.Case(
                json_file,
                error_fixer=error_fixer,
                omim_obj=omim_obj,
                exclude_benign_variants=exclusion
            )

    case_objs = yield_cases(
        jsons,
        error_fixer,
        omim_obj,
        config_data.preprocess.getboolean("exclude_normal_variants")
    )

    mutalyzer.correct_reference_transcripts(case_objs)

    if config_data.general.getboolean('dump_intermediate') and not args.single:
         pickle.dump(case_objs, open('case_cleaned.p', 'wb'))

    return case_objs


def phenomize(args, config_data, cases):
    '''Phenomization using charite phenomization service.'''
    print("== Phenomization of cases ==")
    if "phenomizer" in config_data and config_data.phenomizer["url"]:

        @progress_bar("Phenomization")
        def yield_phenomized(case_objs, phen):
            '''Phenomize a single case. Modifications are inplace.'''
            for case_obj in case_objs:
                case_obj.phenomize(phen)
                yield

        phen = phenomizer.PhenomizerService(config=config_data)
        yield_phenomized(cases, phen)

        if config_data.general.getboolean('dump_intermediate') \
                and not args.single:
            pickle.dump(cases, open('case_phenomized.p', 'wb'))
    else:
        print("No config found. Phenomization will be skipped.")

    return cases


def convert_to_old_format(args, config_data, cases):
    '''Convert case files to old json format objects.'''
    print("== Mapping to old json format ==")
    destination = args.output or config_data.conversion["output_path"]

    omim_obj = omim.Omim(config=config_data)

    @progress_bar("Convert old")
    def yield_old_json(case_objs, destination, omim_obj):
        '''Create an old case object.'''
        for case_obj in case_objs:
            old = json_parser.OldJson.from_case_object(
                case_obj,
                destination,
                omim_obj
            )
            old.save_json()
            yield old

    return yield_old_json(cases, destination, omim_obj)


def get_qc_cases(config_data, cases):
    '''Get qc results for all cases.'''
    omim_obj = omim.Omim(config=config_data)
    return {c.case_id: (c.check(omim_obj), c) for c in cases}


def save_vcfs(args, config_data, qc_cases):
    '''Create VCF files from genetic information and create a config.yml
    listing all vcf files.
    '''
    simulated = config_data.vcf["simulated"]
    realvcf = config_data.vcf["realvcf"]
    config_path = config_data.vcf["config_file"]

    @progress_bar("Generate VCFs")
    def yield_vcf(case_objs):
        '''Dump simulated vcf files.'''
        for case_obj in case_objs:
            case_obj.dump_vcf(simulated)
            yield

    yield_vcf([v[1] for v in qc_cases.values() if v[0][0]])

    if config_data.general.getboolean('dump_intermediate') and not args.single:
        pickle.dump(qc_cases, open('qc_case_with_simulated_vcf.p', 'wb'))

    # only merge if updating changed and processing single cases
    merge = config_data.general.getboolean("changed_only") or args.single
    create_config(config_path, simulated, realvcf, merge)
    return qc_cases


def quality_check_cases(args, config_data, qc_cases, old_jsons):
    '''Output quality check summaries.'''
    print("== Quality check ==")

    # Cases failing qc altogether
    qc_failed_msg = {
        v[1].case_id: v[0] for v in qc_cases.values() if not v[0][0]
    }
    # Cases passing quality check
    qc_passed = {v[1].case_id: v[1] for v in qc_cases.values() if v[0][0]}

    qc_vcf = {
        case_id: case.check_vcf() for case_id, case in qc_passed.items()
    }

    qc_vcf_failed = {
        case_id: result for case_id, result in qc_vcf.items() if not result[0]
    }

    # cases have to pass vcf check
    qc_passed = {
        case_id: case
        for case_id, case in qc_passed.items() if qc_vcf[case_id][0]
    }

    # Cases with mutations marked as benign excluded from analysis
    qc_benign_passed = {
        k: v for k, v in
        {k: v.get_benign_excluded() for k, v in qc_passed.items()}.items()
        if v > 0
    }

    omim_obj = omim.Omim(config=config_data)

    # Cases where pathogenic diagnosed mutation is not in geneList
    @progress_bar("Get pathogenic genes in geneList")
    def pathogenic_genes_process(cases):
        '''Get boolean value, whether pathogenic gene is contained in
        genes converted from detected syndromes.'''
        for case_id, case_obj in cases.items():
            yield case_id, case_obj.pathogenic_gene_in_gene_list(omim_obj)

    qc_pathongenic_passed = dict(
        c for c in pathogenic_genes_process(qc_passed) if not c[1][0]
    )

    # Compiled stats to be dumped into a json file
    qc_output = {
        "failed": qc_failed_msg,
        "benign_excluded": qc_benign_passed,
        "pathogenic_missing": qc_pathongenic_passed,
        "vcf_failed": qc_vcf_failed,
        "passed": {k: "" for k in qc_passed}
    }

    qc_detailed_path = config_data.quality["qc_detailed_log"]
    if config_data.general.getboolean("changed_only") \
            and os.path.exists(qc_detailed_path):
        with open(qc_detailed_path, "r") as qc_file:
            olddata = json.load(qc_file)
        qc_output = {
            k: {**olddata[k], **v}
            for k, v in qc_output.items()
        }

    # save qc results in detailed log if needed
    print("Saving qc log")
    if config_data.quality.getboolean("qc_detailed") \
            and config_data.quality["qc_detailed_log"]:
        with open(qc_detailed_path, "w") as qc_out:
            json.dump(qc_output, qc_out, indent=4)

    # move cases to qc directory
    print("Saving passing cases to new location")
    if config_data.quality["qc_output_path"] and old_jsons:
        # create output directory if needed
        os.makedirs(config_data.quality["qc_output_path"], exist_ok=True)

        old_jsons = {j.get_case_id(): j for j in old_jsons}

        @progress_bar("Save passing qc")
        def save_old_to_qc():
            '''Save old jsons passing QC to a new location.'''
            for pcase in qc_passed:
                old_js = old_jsons[pcase]
                old_js.save_json(
                    destination=config_data.quality["qc_output_path"]
                )
        save_old_to_qc()

    return {"pass": len(qc_passed), "fail": len(qc_failed_msg)}, qc_passed


def main():
    '''
    Some program blocks are enabled and disabled via config options in general
    '''

    configure_logging("lib")
    config_data = config.ConfigManager()

    args = parse_arguments()
    if not args.pickle:
        jsons = create_jsons(args, config_data)
        cases = create_cases(args, config_data, jsons)
    else:
        with open(args.pickle, "rb") as pickled_file:
            cases = pickle.load(pickled_file)

    if args.entry == "pheno":
        cases = phenomize(args, config_data, cases)

    if args.entry == "pheno" or args.entry == "convert":
        old_jsons = convert_to_old_format(args, config_data, cases)
    else:
        old_jsons = None

    if args.entry != "qc":
        # QC Check for only using cases passing qc
        qc_cases = get_qc_cases(config_data, cases)

        if not args.skip_vcf:
            # VCF Generation
            qc_cases = save_vcfs(args, config_data, qc_cases)
    else:
        qc_cases = cases

    # Quality check
    stats, qc_cases = quality_check_cases(
        args, config_data, qc_cases, old_jsons
    )
    print(
        "== QC results ==\nPassed: {pass} Failed: {fail}".format(
            **stats)
    )


if __name__ == '__main__':
    main()
