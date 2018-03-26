# standard libraries
import os
import logging
import pickle

# 3rd party libraries
# import hgvs

# own libraries
from lib import download
from lib.api.phenomizer import PhenomizerService
from lib.api.omim import Omim
from lib.api.mutalyzer import correct_reference_transcripts

from lib.model.json import NewJson, OldJson
from lib.model.case import Case
from lib.model.config import ConfigManager

import hgvs.dataproviders.uta


def main():
    '''Create json files in the old format from new json files.
    '''
    logfile_path = 'logs.txt'
    logging.basicConfig(filename=logfile_path, level=logging.WARNING)
    config = ConfigManager()
    # Load configuration and initialize API bindings
    omim = Omim(config=config)
    phen = PhenomizerService(config=config)

    # Download new files from AWS Bucket
    download.backup_s3_folder(config=config)

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

    print('Unfiltered', len(new_json_objs))

    filtered_new = [j for j in new_json_objs if j.check()[0]]

    print('Filtered rough criteria', len(filtered_new))

    case_objs = [Case(j) for j in filtered_new]
    case_objs = [c for c in case_objs if c.variants]
    print('Cases with created hgvs objects', len(case_objs))

    # batch call mutalyzer
    correct_reference_transcripts(case_objs)
    pickle.dump(case_objs, open('case_cleaned.p', 'wb'))

    #dump vcfs
    #case_objs = pickle.load(open('case_cleaned.p', 'rb'))
    for c in case_objs:
        c.dump_vcf("data/PEDIA/mutations/")


if __name__ == '__main__':
    main()
