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

# own libraries
from lib import download, errorfixer
from lib.visual import progress_bar
from lib.model import json, case, config
from lib.api import phenomizer, omim, mutalyzer


def configure_logging(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    stdout_channel = logging.StreamHandler()
    stdout_channel.setLevel(logging.INFO)
    formatter = logging.Formatter("%(message)s")
    stdout_channel.setFormatter(formatter)
    logger.addHandler(stdout_channel)


def parse_arguments():
    parser = ArgumentParser(description=(
        "Process f2g provided jsons into a format processable by "
        "classification."))
    parser.add_argument("-s", "--single", help="Process a single json file.")
    parser.add_argument("-o", "--output",
                        help="Destination of created old json.",
                        default="")
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

    return json_files, corrected, "process/convert"


@progress_bar("Process jsons")
def yield_jsons(json_files, corrected):
    for json_file in json_files:
        yield json.NewJson.from_file(json_file, corrected)


@progress_bar("Create cases")
def yield_cases(json_files, error_fixer):
    for json_file in json_files:
        yield case.Case(json_file, error_fixer=error_fixer)


@progress_bar("Phenomization")
def yield_phenomized(case_objs, phen):
    for case_obj in case_objs:
        case_obj.phenomize(phen)
        yield


@progress_bar("Convert old")
def yield_old_json(case_objs, destination, omim_obj):
    for case_obj in case_objs:
        old = json.OldJson.from_case_object(case_obj, destination, omim_obj)
        old.save_json()
        yield


def main():
    '''
    Some program blocks are enabled and disabled via config options in general
    '''

    configure_logging("lib")
    config_data = config.ConfigManager()

    # Load configuration and initialize API bindings
    phen = phenomizer.PhenomizerService(config=config_data)
    error_fixer = errorfixer.ErrorFixer(config=config_data)

    args = parse_arguments()

    # get either from single file or from directory
    json_files, corrected, destination = ([args.single], "", args.output) \
        if args.single else json_from_directory(config_data)
    new_json_objs = yield_jsons(json_files, corrected)
    # [json.NewJson.from_file(f, corrected) for f in json_files]
    # for js in new_json_objs:
    #     print(js.get_case_id())
    #     print(js.get_vcf())
    # return

    print('Unfiltered', len(new_json_objs))

    filtered_new = [j for j in new_json_objs if j.check()[0]]

    print('Filtered rough criteria', len(filtered_new))

    case_objs = yield_cases(filtered_new, error_fixer)
    # [case.Case(j, error_fixer=error_fixer) for j in filtered_new]
    case_objs = [c for c in case_objs if c.variants]
    print('Cases with created hgvs objects', len(case_objs))

    mutalyzer.correct_reference_transcripts(case_objs)

    if config_data.general['dump_intermediate']:
        pickle.dump(case_objs, open('case_cleaned.p', 'wb'))

    yield_phenomized(case_objs, phen)
    # for case_obj in case_objs:
    #     case_obj.phenomize(phen)

    if config_data.general['dump_intermediate']:
        pickle.dump(case_objs, open('case_phenomized.p', 'wb'))

    omim_obj = omim.Omim(config=config_data)
    yield_old_json(case_objs, destination, omim_obj)
    # for case_obj in case_objs:
    #     print('Converting to old json', case_obj.case_id)
    #     old_js = json.OldJson.from_case_object(case_obj, destination,
    #                                            omim_obj)
    #     old_js.save_json()


if __name__ == '__main__':
    main()
