# standard libraries
import os
import logging
from pprint import pprint
import json

import pickle

# 3rd party libraries
import hgvs

# own libraries
from lib import download
from lib.api.phenomizer import PhenomizerService
from lib.api.omim import Omim

from lib.model.json import NewJson, OldJson
from lib.model.case import Case
from lib.model.config import ConfigManager

def main():
    logfile_path = 'logs.txt'
    logging.basicConfig(filename=logfile_path,level=logging.WARNING)
    config = ConfigManager()

    ## Load configuration and initialize API bindings
    aws_param = config.preprocess['aws_access_key', 'aws_secret_key', 'download_location']
    omim_param = config.phenomization['OMIM_Key','OMIM_Dir','OMIM_Cached']
    pheno_args = {
            'url':config.phenomization['Phenomizer_Url'],
            'user':config.phenomization['Phenomizer_User'],
            'password':config.phenomization['Phenomizer_Password']
            }
    omim = Omim(*omim_param.values())
    annot = PhenomizerService(**pheno_args)

    ### Download new files from AWS Bucket
    # download.backup_s3_folder(**aws_param)

    ### Initial Quality check of new json
    unprocessed_jsons = os.path.join(aws_param['download_location'],'cases')
    json_files = [os.path.join(unprocessed_jsons,x) for x in os.listdir(unprocessed_jsons) if os.path.splitext(x)[1] == '.json']

    # corrected is a directory which can contain manually edited case jsons
    # that should differ from the original only in content, not in overall structure.
    # this should make resolving some exotic errors a lot easier
    corrected = config.preprocess['corrected_location']
    new_json_objs = [NewJson(f, corrected) for f in json_files]

    print('Unfiltered', len(new_json_objs))

    filtered_new = [ j for j in new_json_objs if j.check()[0] ]
    print('Filtered rough criteria', len(filtered_new))

    case_objs = [ Case(j) for j in filtered_new ]
    case_objs = [ c for c in case_objs if c.variants ]
    print('Cases with created hgvs objects', len(case_objs))

    [ e.syndromes_to_genes(omim) for e in case_objs ]

    pickle.dump(case_objs, open('case_cleaned.p', 'wb'))

    # output_folder = 'curated'
    # os.makedirs(output_folder,exist_ok=True)

    # hdp=hgvs.dataproviders.uta.connect()
    # newQC.checkjson(os.path.join(unprocessed_jsons,json_files[0]), hdp)

if __name__ == '__main__':
    main()
