'''
Bindings for Mutalyzer. A service which can be used to look up RS numbers and
check the validity of hgvs strings.
'''
import logging
import time
import io
import re
from typing import Union, List
import pandas
import zeep
import hgvs.parser
# import pandas

LOGGER = logging.getLogger(__name__)

# alternative hgvs transcript
RE_VERSION_ALTERNATIVE = re.compile('We found these versions: ([\w.]+)')


def correct_reference_transcripts(case_objs: List['Case']) -> List['Case']:
    '''Check reference transcript number via batch call to mutalyzer.
    This will edit the hgvs objects in-place.
    '''
    case_dict = {v.case_id: v.variants for v in case_objs}
    mutalyzer = Mutalyzer()
    case_dict = mutalyzer.correct_transcripts(case_dict)


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
        super().__init__(self.wsdl_url)

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
        LOGGER.info("Submitting batch job to mutalyzer.")
        while True:
            remaining_jobs = self.service.monitorBatchJob(batch_id)
            LOGGER.info('Remaining %d', remaining_jobs)
            if remaining_jobs == 0:
                break
            time.sleep(1)
        LOGGER.info("Finished batch job.")
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
        # create transcript input data
        data = "\n".join([str(v) for l in transcripts.values() for v in l])
        response = self.batch_position_convert(data=data)
        for hgvs_variants in transcripts.values():
            for var in hgvs_variants:
                key = str(var)
                if key in response:
                    alt_transcript = response[key]
                    if alt_transcript:
                        LOGGER.info('Replace %s with %s',
                                    var.ac, alt_transcript)
                        var.ac = alt_transcript
        return transcripts

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
