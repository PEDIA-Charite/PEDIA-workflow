'''
Bindings for Mutalyzer. A service which can be used to look up RS numbers and
check the validity of hgvs strings.
'''
import logging
import time
import io
import os
import json
import re
from typing import Union, List
import pandas
import requests
from requests.adapters import HTTPAdapter
import zeep
from zeep.transports import Transport

from lib import visual
from lib.constants import CACHE_DIR

LOGGER = logging.getLogger(__name__)

# alternative hgvs transcript
RE_VERSION_ALTERNATIVE = re.compile(r'We found these versions: ([\w.]+)')


def correct_reference_transcripts(case_objs: List['Case']) -> List['Case']:
    '''Check reference transcript number via batch call to mutalyzer.
    This will edit the hgvs objects in-place.
    '''
    case_dict = {
        v.case_id: [
            vv for m in v.get_hgvs_models() if not m.corrected
            for vv in m.variants
        ]
        for v in case_objs
    }
    mutalyzer = Mutalyzer()
    mutalyzer.correct_transcripts(case_dict)


def check_errors(errordata) -> Union[str, None]:
    '''Check whether transcript number contains errors and try to fix these.
    '''
    if pandas.notna(errordata):
        match = RE_VERSION_ALTERNATIVE.search(errordata)
        if match:
            return match.group(1)
    return None


class Mutalyzer(zeep.Client):
    '''Implements API bindings for the Mutalyzer.
    '''

    base_url = 'https://mutalyzer.nl/json/'
    wsdl_url = 'https://mutalyzer.nl/services/?wsdl'

    def __init__(self):
        session = requests.Session()
        session.mount('http', HTTPAdapter(max_retries=3))
        session.mount('https', HTTPAdapter(max_retries=3))
        transport = Transport(session=session)

        super().__init__(self.wsdl_url, transport=transport)

        self._transcript_cache = self._load_cache()

    def batch_position_convert(self, data: str):
        '''Submit a batch job to the mutalyzer, monitor it and return the
        output data.
        '''
        # data needs to be base64 encoded
        # data = base64.b64encode(data.encode('utf-8'))
        data = data.encode('utf-8')
        method = 'PositionConverter'
        # used to denote genome build when using PositionConverter
        # alternatively use GRCh37
        argument = 'GRCh37'

        batch_id = self.service.submitBatchJob(
            data=data, process=method, argument=argument)

        # wait for the batch job to finish
        LOGGER.debug("Submitting batch job to mutalyzer.")
        max_obj = 0
        cur_obj = 0
        while True:
            remaining_jobs = self.service.monitorBatchJob(batch_id)
            max_obj = max(remaining_jobs, max_obj)
            cur_obj = max_obj - remaining_jobs
            LOGGER.debug('Remaining %d', remaining_jobs)
            visual.print_status("Mutalyzer", width=20,
                                cur=cur_obj, size=max_obj)
            if remaining_jobs == 0:
                print("")
                break
            time.sleep(1)
        LOGGER.debug("Finished batch job.")
        # get the batch job results
        batch_result = self.service.getBatchJob(batch_id)
        result_string = batch_result.decode('utf-8')
        # returned tsv file is really weird
        # first three columns are single entry, while the last col gets
        # all tsv cells until the next row
        names = ['Input', 'Errors', 'Chromosomal', 'Codings']
        result_table = pandas.read_table(
            io.StringIO(result_string), sep='\t',
            names=names, index_col=0)
        result_table.drop(index='Input Variant', inplace=True)
        transcript_alt = result_table['Errors'].apply(check_errors)

        return transcript_alt.to_dict()

    def correct_transcripts(self, transcripts: dict) -> dict:
        '''Get a dictionary with assignments for each transcript we have sent
        to the batch job, whether the given transcript is correct.
        '''
        # search in cache first
        all_variants = [v for l in transcripts.values() for v in l]
        remaining_transcripts = [
            v for v in all_variants if not self._modify_transcript_cached(v)
        ]

        # create transcript input data
        data = "\n".join([str(v) for v in remaining_transcripts])
        if not data:
            LOGGER.warning("Data empty. No batch process created.")
            return []
        response = self.batch_position_convert(data=data)
        for var in remaining_transcripts:
            self._modify_transcript(var, response)
        self._update_cache(response)
        return transcripts

    def _get_cache_path(self) -> str:
        return os.path.join(CACHE_DIR, __name__ + "_transcript_cache.json")

    def _load_cache(self) -> dict:
        if os.path.exists(self._get_cache_path()):
            with open(self._get_cache_path(), "r") as cache_file:
                data = json.load(cache_file)
        else:
            data = {}
        return data

    def _update_cache(self, update: dict) -> None:
        self._transcript_cache = {**self._transcript_cache, **update}
        os.makedirs(CACHE_DIR, exist_ok=True)
        with open(self._get_cache_path(), "w") as cache_file:
            json.dump(self._transcript_cache, cache_file)

    def _modify_transcript_cached(
            self,
            variant: "hgvs.sequencevariant"
    ) -> bool:
        return self._modify_transcript(variant, self._transcript_cache)

    def _modify_transcript(
            self,
            variant: "hgvs.sequencevariant",
            trans_dict: dict
    ) -> bool:
        if str(variant) in trans_dict:
            alt_transcript = trans_dict[str(variant)]
            if alt_transcript:
                LOGGER.debug(
                    'Replace %s with %s', variant.ac, alt_transcript
                )
                variant.ac = alt_transcript
            return True
        return False

    def get_db_snp_descriptions(self, rs_id: str) -> [str]:
        '''Return a list of possible RS numbers for the given RS code.
        '''
        tries = 0
        while tries < 5:
            try:
                variants = self.service.getdbSNPDescriptions(rs_id)
                break
            except zeep.exceptions.Fault as error:
                print(error)
                tries += 1
        return variants

    def check_syntax(self, hgvs_string: str) -> dict:
        '''Check the syntax of an hgvs variant and return the validity and list
        of possible errors.
        '''
        return self.service.checkSyntax(hgvs_string)
