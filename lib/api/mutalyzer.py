'''
Bindings for Mutalyzer. A service which can be used to look up RS numbers and
check the validity of hgvs strings.
'''
import time
import io
import re
from typing import Union, List
import pandas
import zeep
import hgvs.parser
# import pandas

# alternative hgvs transcript
RE_VERSION_ALTERNATIVE = re.compile('We found these versions: ([\w.]+)')

HGVS_PARSER = hgvs.parser.Parser()


def correct_reference_transcripts(case_objs: List['Case']) -> List['Case']:
    '''Check reference transcript number via batch call to mutalyzer.
    This will edit the hgvs objects in-place.
    '''
    case_dict = {v.case_id: v.variants for v in case_objs}
    mutalyzer = Mutalyzer()
    case_dict = mutalyzer.correctTranscripts(case_dict)


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

    def batchPositionConvert(self, data: str):
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
        while True:
            remaining_jobs = self.service.monitorBatchJob(batch_id)
            print('Remaining', remaining_jobs)
            if remaining_jobs == 0:
                break
            time.sleep(1)
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

    def correctTranscripts(self, transcripts: dict) -> dict:
        '''Get a dictionary with assignments for each transcript we have sent
        to the batch job, whether the given transcript is correct.
        '''
        # create transcript input data
        data = "\n".join([str(v) for l in transcripts.values() for v in l])
        response = self.batchPositionConvert(data=data)
        for entry_id, hgvs_variants in transcripts.items():
            for var in hgvs_variants:
                key = str(var)
                if key in response:
                    alt_transcript = response[key]
                    if alt_transcript:
                        print('Replace {} with {}'.format(
                            var.ac, alt_transcript))
                        var.ac = alt_transcript
        return transcripts

    def getdbSNPDescriptions(self, rs_id: str) -> [str]:
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

    def checkSyntax(self, hgvs_string: str) -> dict:
        '''Check the syntax of an hgvs variant and return the validity and list
        of possible errors.
        '''
        return self.service.checkSyntax(hgvs_string)


def main():
    '''Only for testing purposes.
    '''
    mut = Mutalyzer()
    print('Convert RS number to hgvs variant')
    hgvs_list = mut.getdbSNPDescriptions("rs386834107")
    print(hgvs_list)
    print('Check HGVS syntax. This does not check transcript versions')
    check = mut.checkSyntax('NM_001127178.2:c.2005C>T')
    print(check)

    print('Regenerate missing transcript version number')
    hgvs_strings = {
        'yolo': 'NM_003002.3:c.274G>T',
        'foker': 'LRG_9t1:c.274G>T',
        'poker': 'chr11:g.111959693G>T',
        'condor': 'NC_000011.9:g.111959693G>T',
        'missing_version': 'NM_001127178:c.2005C>T'
    }

    hgvs_strings = {k: [HGVS_PARSER.parse_hgvs_variant(v)]
                    for k, v in hgvs_strings.items()}
    mut.correctTranscripts(hgvs_strings)
    print(hgvs_strings)


if __name__ == '__main__':
    main()
