'''
Convert new jsons to old jsons
---
This script has been used to convert information in the new jSON format to the
older json format.

This script should only be used as a reference for implementation and could be
outdated for later API changes.
'''
# standard libraries
import os
import logging

import pickle
import sys

sys.path.append(os.getcwd())

# 3rd party libraries
# import hgvs

# own libraries
# from lib import download
from lib.api.phenomizer import PhenomizerService
from lib.api.omim import Omim
from lib.api.mutalyzer import correct_reference_transcripts

from lib.model.json import NewJson, OldJson
from lib.model.case import Case
from lib.model.config import ConfigManager


def main():
    '''Create json files in the old format from new json files.
    '''
    logfile_path = 'logs.txt'
    logging.basicConfig(filename=logfile_path, level=logging.WARNING)
    config = ConfigManager()

    # Load configuration and initialize API bindings
    # omim = Omim(config=config)
    # phen = PhenomizerService(config=config)

    # Download new files from AWS Bucket
    # download.backup_s3_folder(config=config)

    # Initial Quality check of new json
    unprocessed_jsons = os.path.join(config.aws['download_location'], 'cases')
    json_files = [os.path.join(unprocessed_jsons, x)
                  for x in os.listdir(unprocessed_jsons)
                  if os.path.splitext(x)[1] == '.json']

    # corrected is a directory which can contain manually edited case jsons
    # that should differ from the original only in content, not in overall
    # structure.  this should make resolving some exotic errors a lot easier
    corrected = config.preprocess['corrected_location']
    new_json_objs = [NewJson.from_file(f, corrected) for f in json_files]

    # print('Unfiltered', len(new_json_objs))
    filtered_new = [j for j in new_json_objs if j.check()[0]]
    vcf_directory = os.path.join(config.aws['download_location'], "vcfs")
    for json_case in filtered_new:
        json_case.get_vcf(vcf_directory)


if __name__ == '__main__':
    main()
