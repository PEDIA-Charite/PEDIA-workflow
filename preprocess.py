#!/usr/bin/env python3
'''
Preprocessing script running the preprocessing based on configuration
options.
'''
# standard libraries
import os
import logging
import logging.config

import pickle

# own libraries
from lib import download, errorfixer
from lib.model import json, case, config
from lib.api import phenomizer, omim, mutalyzer


def main():
    '''
    Some program blocks are enabled and disabled via config options in general
    '''
    # configure logging
    # logfile_path = 'logs.txt'
    logger = logging.getLogger("lib")
    logger.setLevel(logging.INFO)
    stdout_channel = logging.StreamHandler()
    stdout_channel.setLevel(logging.INFO)
    formatter = logging.Formatter("%(message)s")
    stdout_channel.setFormatter(formatter)
    logger.addHandler(stdout_channel)

    config_data = config.ConfigManager()

    # Load configuration and initialize API bindings
    phen = phenomizer.PhenomizerService(config=config_data)
    error_fixer = errorfixer.ErrorFixer(config=config_data)

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
    new_json_objs = [json.NewJson.from_file(f, corrected) for f in json_files]

    print('Unfiltered', len(new_json_objs))

    filtered_new = [j for j in new_json_objs if j.check()[0]]

    print('Filtered rough criteria', len(filtered_new))

    case_objs = [case.Case(j, error_fixer=error_fixer) for j in filtered_new]
    case_objs = [c for c in case_objs if c.variants]
    print('Cases with created hgvs objects', len(case_objs))

    mutalyzer.correct_reference_transcripts(case_objs)

    if config_data.general['dump_intermediate']:
        pickle.dump(case_objs, open('case_cleaned.p', 'wb'))

    for case_obj in case_objs:
        print('Phenomizing', case_obj.case_id)
        case_obj.phenomize(phen)

    if config_data.general['dump_intermediate']:
        pickle.dump(case_objs, open('case_phenomized.p', 'wb'))

    omim_obj = omim.Omim(config=config_data)
    for case_obj in case_objs:
        print('Converting to old json', case_obj.case_id)
        old_js = json.OldJson.from_case_object(case_obj, 'process/convert',
                                               omim_obj)
        old_js.save_json()


if __name__ == '__main__':
    main()
