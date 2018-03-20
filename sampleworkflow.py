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

    #check for format violations
    filtered_new = [j for j in new_json_objs if j.check()[0]]

    #intilize dataprovider for variant mapping
    hdp = hgvs.dataproviders.uta.connect()

    #create case objects
    case_objs = [Case(j,hdp) for j in filtered_new]

    #number of cases that contain variants
    #case_objs = [c for c in case_objs if c.variants]
    #print(len(case_objs))

    #pickle.dump(case_objs, open('case_cleaned.p', 'wb'))
    #case_objs = pickle.load(open('case_corrected.p', 'rb'))

    #get cases without vcf
    failed_case_objs = [c for c in case_objs if type(c.vcf)==list]
    #try to correct the transcripts using mutalyzer
    correct_reference_transcripts(failed_case_objs)
    #try again to create vcfs
    for c in failed_case_objs:
        c.vcf=c.create_vcf(hdp)

    #pickle.dump(case_objs, open('case_corrected.p', 'wb'))
    #dump vcfs to path VCF
    for c in case_objs:
        c.dump_vcf("./VCF/")

    #dump information about failed HGVS codes to an output file
    #with open("failedhgvs.txt","w+") as output:
    #    output.write("Case\tHGVS Codes\tErrors\n")
    #    for c in case_objs:
    #        if type(c.vcf)==list:
    #            output.write(c.case_id+"\t"+",".join(map(str,c.variants))+"\t"+",".join(map(str,c.vcf))+"\n")


if __name__ == '__main__':
    main()
