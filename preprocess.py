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
import json as json_lib

# own libraries
from lib import download, errorfixer
from lib.visual import progress_bar
from lib.model import json, case, config
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

    return parser.parse_args()


def json_from_directory(config_data: config.ConfigManager) \
        -> Tuple[List[str], str]:
    '''Get a list of json file paths.'''
    # Download new files from AWS Bucket
    if config_data.general["download"]:
        download.backup_s3_folder(config=config_data)

    # Initial Quality check of new json
    unprocessed_jsons = os.path.join(
        config_data.aws['download_location'], 'cases')
    json_files = [os.path.join(unprocessed_jsons, x)
                  for x in os.listdir(unprocessed_jsons)
                  if os.path.splitext(x)[1] == '.json']
    # corrected is a directory which can contain manually edited case jsons
    # that should differ from the original only in content, not in overall
    # structure.
    # this should make resolving some exotic errors a lot easier
    corrected = config_data.preprocess['corrected_location']

    return json_files, corrected

def create_config(simvcffolder: str = "data/PEDIA/mutations", vcffolder: str = "data/PEDIA/vcfs/original") -> None:
    '''Creates config.yml file based on the VCF files'''
    vcffiles = [file.split(".")[0] for file in os.listdir(vcffolder)]
    singlefiles = [file.split(".")[0] for file in os.listdir(simvcffolder)]
    testfiles = []
    with open("config.yml","w") as configfile:
        configfile.write('SINGLE_SAMPLES: \n')
        for file in singlefiles:
            if file not in vcffiles:
                configfile.write(" - " + file + "\n")
        configfile.write('VCF_SAMPLES: \n')
        for file in vcffiles:
            if file in singlefiles:
                configfile.write(" - " + file + "\n")
            else:
                testfiles.append(file)
        configfile.write('TEST_SAMPLES: \n')
        for file in testfiles:
            configfile.write(" - " + file + "\n")


@progress_bar("Process jsons")
def yield_jsons(json_files, corrected):
    for json_file in json_files:
        yield json.NewJson.from_file(json_file, corrected)


@progress_bar("Create cases")
def yield_cases(json_files, error_fixer, exclusion):
    for json_file in json_files:
        yield case.Case(
            json_file,
            error_fixer=error_fixer,
            exclude_benign_variants=exclusion
        )


@progress_bar("Phenomization")
def yield_phenomized(case_objs, phen):
    for case_obj in case_objs:
        case_obj.phenomize(phen)
        yield


@progress_bar("Convert old")
def yield_old_json(case_objs, destination, omim_obj):
    for case_obj in case_objs:
        old = json.OldJson.from_case_object(
            case_obj,
            destination,
            omim_obj
        )
        old.save_json()
        yield old


@progress_bar("Generate VCFs")
def yield_vcf(case_objs, destination):
    for case_obj in case_objs:
        case_obj.dump_vcf(destination)
        yield


def create_jsons(args, config_data):
    print("== Process new json files ==")
    # get either from single file or from directory
    json_files, corrected = ([args.single], "") \
        if args.single else json_from_directory(config_data)
    new_json_objs = yield_jsons(json_files, corrected)

    print('Unfiltered', len(new_json_objs))

    filtered_new = [j for j in new_json_objs if j.check()[0]]
    print('Filtered rough criteria', len(filtered_new))
    return filtered_new


def create_cases(args, config_data, jsons):
    print("== Create cases from new json format ==")
    error_fixer = errorfixer.ErrorFixer(config=config_data)
    case_objs = yield_cases(
        jsons,
        error_fixer,
        config_data.preprocess["exclude_normal_variants"]
    )

    mutalyzer.correct_reference_transcripts(case_objs)

    if config_data.general['dump_intermediate']:
         pickle.dump(case_objs, open('case_cleaned.p', 'wb'))

    return case_objs


def phenomize(config_data, cases):
    print("== Phenomization of cases ==")
    phen = phenomizer.PhenomizerService(config=config_data)
    yield_phenomized(cases, phen)

    if config_data.general['dump_intermediate']:
        pickle.dump(cases, open('case_phenomized.p', 'wb'))
    return cases


def convert_to_old_format(args, config_data, cases):
    print("== Mapping to old json format ==")
    destination = args.output or config_data.conversion["output_path"]

    omim_obj = omim.Omim(config=config_data)
    return yield_old_json(cases, destination, omim_obj)


def save_vcfs(config_data, cases):
    yield_vcf(cases, 'data/PEDIA/mutations')
    cases = [case for case in cases if hasattr(case, 'vcf')]
    if config_data.general['dump_intermediate']:
        pickle.dump(cases, open('case_with_simulated_vcf.p', 'wb'))
    return cases


def quality_check_cases(args, config_data, cases, old_jsons):
    '''Output quality check summaries.'''
    print("== Quality check ==")
    omim_obj = omim.Omim(config=config_data)
    qc_results = {
        c.case_id: c.check(omim_obj)
        for c in cases
    }

    qc_failed = {c: q for c, q in qc_results.items() if not q[0]}

    passed_cases = [
        c for c in cases if c.check(omim_obj)[0]
    ]
    # save qc results in detailed log if needed
    if config_data.quality["qc_detailed"] \
            and config_data.quality["qc_detailed_log"]:
        with open(config_data.quality["qc_detailed_log"], "w") as qc_out:
            json_lib.dump(qc_failed, qc_out, indent=4)

    # move cases to qc directory
    if config_data.quality["qc_output_path"] and old_jsons:
        # create output directory if needed
        os.makedirs(config_data.quality["qc_output_path"], exist_ok=True)

        old_jsons = {j.get_case_id(): j for j in old_jsons}

        @progress_bar("Save passing qc")
        def save_old_to_qc():
            for pcase in passed_cases:
                old_js = old_jsons[pcase.get_case_id()]
                old_js.save_json(
                    destination=config_data.quality["qc_output_path"]
                )


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

    cases = [case for case in cases if case.check()[0]]
    if args.entry == "pheno":
        cases = phenomize(config_data, cases)

    if args.entry == "pheno" or args.entry == "convert":
        old_jsons = convert_to_old_format(args, config_data, cases)
    else:
        old_jsons = None

    if args.entry == "pheno" \
            or args.entry == "convert" \
            or args.entry == "qc":
        quality_check_cases(args, config_data, cases, old_jsons)

    cases = save_vcfs(config_data, cases)
    create_config()


if __name__ == '__main__':
    main()
