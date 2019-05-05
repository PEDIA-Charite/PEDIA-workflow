#!/usr/bin/env python3
'''
Preprocessing script running the preprocessing based on configuration
options.
'''
# standard libraries
import os
import sys
import shutil
import logging
import logging.config
from typing import Tuple, List

import json

import yaml
import snakemake.workflow

# own libraries
from lib import errorfixer, quality_check, pickler
from lib.processor import Processor
from lib.visual import progress_bar, multiprocess
from lib.model import json_parser, case, config, args_parser
from lib.api import omim, mutalyzer, jannovar

from lib.global_singletons import AWS_INST, MUTALYZER_INST, LAB_INST


def configure_logging(logger_name, logger_file: str = "preprocess.log"):
    '''Set up logging devices for logging to screen and a separate file
    with different log levels.'''
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

def create_config(
        config_path: str = "config.yml",
        simvcffolder: str = "data/PEDIA/mutations",
        vcffolder: str = "data/PEDIA/vcfs/original"
) -> None:
    '''Creates config.yml file based on the VCF files'''
    # real vcf files
    vcffiles = {
        int(f.split(".")[0])
        for f in os.listdir(vcffolder) if not f.startswith(".")
    }
    # simulated vcf files
    singlefiles = {
        int(f.split(".")[0])
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
    with open("config.yml", "w") as configfile:
        yaml.dump(config_data, configfile, default_flow_style=False)


def create_jsons(config_data, convert_failed):
    '''Create a list of new formatjson objects.'''
    print("== Process new json files ==")
    # get either from single file or from directory
    corrected = config_data.input["corrected_path"]
    if config_data.input["input_files"]:
        json_files = config_data.input["input_files"]
    elif config_data.input["lab"] and config_data.input["lab_case_id"]:
        json_folder = config_data.input["download_path"]
        unprocessed_jsons = os.path.join(json_folder, 'cases')
        json_files = [LAB_INST.download_lab_case(unprocessed_jsons, config_data.input["lab_case_id"])]
        print("downalod lab case")
    else:
        json_folder = config_data.input["download_path"]
        if config_data.input["download"]:
            LAB_INST.download_all_lab_case(json_folder)

        unprocessed_jsons = os.path.join(json_folder, 'cases')

        json_files = [
            os.path.join(unprocessed_jsons, x)
            for x in os.listdir(unprocessed_jsons)
            if os.path.splitext(x)[1] == '.json'
        ]

    if config_data.input["aws_format"]:
        new_json_objs = progress_bar("Process jsons")(
            lambda x, y: json_parser.NewJson.from_file(x, y)
        )(json_files, corrected)
    else:
        new_json_objs = progress_bar("Process jsons")(
            lambda x, y: json_parser.LabJson.from_file(x, y)
        )(json_files, corrected)

    print('Unfiltered', len(new_json_objs))
    json_failed_data = {
        "json_check_failed": {
            case_id: {
                "issues": issues,
            }
            for case_id, (valid, issues) in
            [(j.get_case_id(), j.check(convert_failed)) for j in new_json_objs]
            if not valid
        }
    }

    filtered_new = [j for j in new_json_objs if j.check(convert_failed)[0]]

    if convert_failed:
        failed_jsons = [j for j in new_json_objs if not j.check(convert_failed)[0]]
    else:
        failed_jsons = []

    print('Filtered rough criteria', len(filtered_new))
    return filtered_new, failed_jsons, json_failed_data


def touch_hgvs(case):
    case.hgvs_models
    return case


def create_cases(config_data, jsons):
    '''Create cases from list of jsons.'''
    print("== Create cases from new json format ==")

    case_objs = progress_bar("Create cases")(
        lambda json_file: case.Case(json_file, config_data, vcf_path=config_data.input["vcf"])
    )(jsons)

    print("Correcting transcripts with mutalyzer")

    case_objs = multiprocess("Fetch hgvs", touch_hgvs, case_objs)

    MUTALYZER_INST.correct_reference_transcripts(case_objs)

    print("Creating pickle.")
    if config_data.dump_intermediate:
        with open('case_cleaned.p', 'wb') as pfile:
            pickler.CasePickler(pfile).dump(case_objs,)

    return case_objs


def create_old_json(case_obj, destination):
    '''Create an old case object.'''
    old = json_parser.OldJson.from_case_object(case_obj, destination)
    old.save_json()
    return old, case_obj


def convert_to_old_format(config_data, cases):
    '''Convert case files to old json format objects.'''
    print("== Mapping to old json format ==")
    result = multiprocess(
        "Create old", create_old_json, cases,
        destination=config_data.output["converted_path"]
    )
    old_jsons = [old for old, _ in result]
    cases = [c for _, c in result]
    return old_jsons, cases

def convert_failed_cases(config_data, jsons):
    '''Convert failed cases to old json format objects'''
    cases = create_cases(config_data, jsons)
    result = multiprocess(
        "Create old", create_old_json, cases,
        destination=config_data.output["converted_path"]
    )

def create_qc_case(case):
    return case.check(), case


def get_qc_cases(config_data, cases):
    '''Get qc results for all cases.'''
    print("== Get QC results for cases ==")
    return multiprocess("QC cases", create_qc_case, cases)


def save_vcfs(config_data, qc_cases):
    '''Create VCF files from genetic information and create a config.yml
    listing all vcf files.
    '''
    simulated = config_data.output["simulated_vcf_path"]
    realvcf = config_data.output["real_vcf_path"]
    config_path = config_data.output["vcf_config_file"]

    progress_bar("Generate VCF")(
        lambda x: x.put_hgvs_vcf(simulated, recreate=False)
    )([case for (valid, _), case in qc_cases if valid])

    print("Pickling vcf cases")
    if config_data.dump_intermediate:
        with open('qc_case_with_simulated_vcf.p', 'wb') as pfile:
            pickler.CasePickler(pfile).dump(qc_cases)

    create_config(config_path, simulated, realvcf)

    return qc_cases


def quality_check_cases(config_data, qc_cases, old_jsons, json_log):
    '''Output quality check summaries.'''
    print("== Quality check ==")

    # Cases failing qc altogether
    qc_failed_msg = {
        case.case_id: (valid, issues)
        for (valid, issues), case in qc_cases if not valid
    }
    # cases with multiple diagnosis
    multi_no_omim = {
        k: v for k, v in
        [
            (case.case_id, case.get_diagnosis()) for _, case in qc_cases
        ]
        if any(str(d["omim_id"]) == "0" for d in v) and len(v) > 1
    }

    # Cases passing quality check
    qc_passed = {case.case_id: case for (valid, _), case in qc_cases if valid}

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

    # Cases where pathogenic diagnosed mutation is not in geneList
    qc_pathongenic_passed = dict(
        (cid, x) for (cid, x) in progress_bar("Pathogenic genes in gene list")(
            lambda c: (c[0], c[1].pathogenic_gene_in_gene_list())
        )(qc_passed.items()) if not x[0]
    )

    # Compiled stats to be dumped into a json file
    qc_output = {
        "failed": qc_failed_msg,
        "benign_excluded": qc_benign_passed,
        "pathogenic_missing": qc_pathongenic_passed,
        "vcf_failed": qc_vcf_failed,
        "multi_no_omim": multi_no_omim,
        "passed": {k: '' for k in qc_passed.keys()},
    }

    qc_output = {**qc_output, **json_log}

    # save qc results in detailed log if needed
    log_path = config_data.output["quality_check_log"]
    if config_data.output["create_log"]:
        # move old file to new location
        if os.path.exists(log_path):
            shutil.move(log_path, log_path+".old")
        print("Saving qc log")
        with open(log_path, "w") as qc_out:
            json.dump(qc_output, qc_out, indent=4)

    # move cases to qc directory
    if old_jsons:
        print("Saving passing cases to new location")
        # create output directory if needed
        os.makedirs(config_data.output["valid_case_path"], exist_ok=True)

        old_jsons = {j.get_case_id(): j for j in old_jsons}

        @progress_bar("Save passing qc")
        def save_old_to_qc(case_obj):
            '''Save old jsons passing QC to a new location.'''
            old_js = old_jsons[case_obj]
            old_js.save_json(save_path=config_data.output["valid_case_path"])

        save_old_to_qc(qc_passed)

    if not config_data.output["create_log"]:
        print(json.dumps(qc_output, indent=4))

    return {
        "pass": len(qc_passed),
        "fail": len(qc_failed_msg) + len(qc_vcf_failed)
    }, qc_passed

def run_workflow(case_id, config_data):
    print("== Start PEDIA workflow == ")
    snakefile = 'Snakefile'
    target_file = os.path.join(config_data.output['output_path'], 'results', str(case_id), 'run.out')
    snakemake.snakemake(snakefile, targets=[target_file], workdir='.')
    print("== PEDIA workflow is completed == ")

def main():
    '''
    Some program blocks are enabled and disabled via config options in general
    '''

    args = args_parser.PEDIAParser().args

    config_data = config.PEDIAConfig(args)
    configure_logging("lib", config_data.logfile_path)

    json_log = {}

    if not args.pickle:
        jsons, failed_jsons, json_log = create_jsons(config_data, args.convert_failed)

        with open("failed_cases.json","w") as failed_log:
            json.dump(json_log,failed_log)

        cases = create_cases(config_data, jsons)
    else:
        with open(args.pickle, "rb") as pickled_file:
            cases = pickler.CaseUnpickler(pickled_file).load()

    if args.convert_failed:
        print("== Convert failed cases == ")
        convert_failed_cases(config_data,failed_jsons)

    if args.entry == "pheno" or args.entry == "convert":
        old_jsons, cases = convert_to_old_format(config_data, cases)
    else:
        old_jsons = None

    if args.entry != "qc":
        # QC Check for only using cases passing qc
        qc_cases = get_qc_cases(config_data, cases)

        if not args.skip_vcf:
            # VCF Generation
            qc_cases = save_vcfs(config_data, qc_cases)
    else:
        qc_cases = cases

    # Quality check
    stats, qc_cases = quality_check_cases(
        config_data, qc_cases, old_jsons, json_log
    )

    print(
        "== QC results ==\nPassed: {pass} Failed: {fail}".format(
            **stats)
    )

    if not args.single and not args.lab:
        quality_check.diff_quality_check(
            config_data.output["quality_check_log"]
        )

    if args.vcf:
        run_workflow(cases[0].case_id, config_data)

if __name__ == '__main__':
    main()
