'''
Phenomization
---
Using a list of HPO features a list of associated diseases with causal genes
and a likelyhood score is returned.

Because our current approach is on a gene-level scoring. This information will
be rearranged into a gene-to-score mapping. Thus the associated diseases
returned by phenomization and boqa can be different.

'''
import io
import re
import os
import logging

import requests
import requests_cache
import pandas

from lib.constants import CACHE_DIR


RE_SYMBOL = re.compile('(\w+) \(\d+\)')
LOGGER = logging.getLogger(__name__)


def match_symbol(descriptor: str) -> str:
    '''Extract the first string from the disease name or gene symbol field.
    Both are formatted so that the respective id is included with the string
    in the format: NAME (ID)
    '''
    match = RE_SYMBOL.search(descriptor)
    # try to match the pattern, if no pattern has been found fall back to
    # returning the original string
    return match and match.group(1) or descriptor


def get_max_gene(df_group: pandas.DataFrame) -> pandas.Series:
    '''Get the gene entry with the highest value.
    '''
    i = df_group['value'].idxmax()
    row = df_group.loc[i]
    return row


class PhenomizerService(requests_cache.CachedSession):
    '''Handling of interop with Phenomizer service, which provides the pheno
    and boqa scores used in the process.
    '''

    def __init__(
            self,
            Phenomizer_Url: str = '',
            Phenomizer_User: str = '',
            Phenomizer_Password: str = '',
            config: 'ConfigParser' = None
    ):
        '''
        Create a new phenomizer service instance.

        Params:
            Phenomizer_Url: Url of phenomizer service
            Phenomizer_User: Username for the service
            Phenomizer_Password: Password for the service
            config: Alternative ConfigParser object to fill url, user and
                    password, which will read the values from a config.ini
        '''
        os.makedirs(CACHE_DIR, exist_ok=True)
        cache_name = os.path.join(CACHE_DIR, __name__+"cache.sqlite")
        super().__init__(cache_name=cache_name, expire_after=None)
        if config:
            self.url = config.phenomizer['url']
            self.user = config.phenomizer['user']
            self.password = config.phenomizer['password']
        else:
            self.url = Phenomizer_Url
            self.user = Phenomizer_User
            self.password = Phenomizer_Password

        # retries settings to repeat api calls in case of failure
        retry = requests.packages.urllib3.util.retry.Retry(
            total=3, read=3, connect=3, backoff_factor=0.3,
            status_forcelist=(500, 404))
        adapter = requests.adapters.HTTPAdapter(max_retries=retry)
        self.mount('http://', adapter)
        self.mount('https://', adapter)

    def _request_data_as_df(self, params: dict, names: list,
                            prefilter: {str: str}) -> pandas.DataFrame:
        '''Get information using api calls and convert to pandas dataframe.
        Args:
            prefilter: Set whether dataframes should be filtered down to
                       disease-id containing OMIM.
        '''
        response = self.get(self.url, params=params)
        # explicit error for faulty request
        response.raise_for_status()
        # remove lines starting with # or empty lines from the returned text
        # response
        rstring = "\n".join(
            [s for s in response.text.split('\n')
             if not s.startswith('#') or s == ''])
        rawdata = io.StringIO(rstring)
        dataframe = pandas.read_table(
            rawdata, sep='\t', index_col=None, header=None, names=names)
        if prefilter:
            for colname, required_string in prefilter.items():
                dataframe = dataframe.loc[
                    dataframe[colname].str.contains(required_string)]
        return dataframe

    def disease_boqa_phenomize(self, hpo_ids: [str]) -> pandas.DataFrame:
        '''Get phenomizer and boqa scorings for the list of hpo ids. A datafame
        joined on the syndrome omim id will be returned.
        '''
        if not hpo_ids:
            scaffold = {
                'disease-id_pheno': "int",
                'disease-id_boqa': "int",
                'value_pheno': "float",
                'value_boqa': "float",
                'gene-id': "object",
                'gene-symbol': "object"
            }
            empty = pandas.DataFrame(columns=scaffold.keys())
            empty = empty.astype(dtype=scaffold)
            return empty
        # this process might need to be retried, but it is currently reliable
        # enough to run directly
        hpo_df = self._request_phenomize(
            hpo_ids, prefilter={'disease-id': 'OMIM'})
        boqa_df = self._request_boqa(
            hpo_ids, prefilter={'disease-id': 'OMIM'})

        # preprocess result for join based on omim disease id
        hpo_df['value'] = 1 - hpo_df['value']
        # remove id tags
        hpo_df['disease-id'] = hpo_df['disease-id'].apply(
            lambda x: x.split(':')[-1]).astype(str)
        hpo_df = hpo_df.set_index('disease-id')
        boqa_df['disease-id'] = boqa_df['disease-id'].apply(
            lambda x: x.split(':')[-1]).astype(str)
        boqa_df = boqa_df.set_index('disease-id')
        # join dataframes on disease id
        scores_df = hpo_df.join(
            boqa_df, how='outer', lsuffix='_pheno', rsuffix='_boqa')
        return scores_df

    def _request_phenomize(self, hpo_ids: [str], prefilter: {str: str}={}) \
            -> pandas.DataFrame:
        '''Get phenomizer information on a list of omim ids.
        '''
        hpo_ids = ','.join(hpo_ids)
        phen_names = ['value', 'score', 'disease-id',
                      'disease-name', 'gene-symbol', 'gene-id']
        params = {
            'mobilequery': 'true',
            'username': self.user,
            'password': self.password,
            'terms': hpo_ids,
            'numres': 100
                }
        dataframe = self._request_data_as_df(
            params=params, names=phen_names, prefilter=prefilter)
        # the score column is not used in the phenomizer for our usecase
        dataframe = dataframe.drop(['score'], axis=1)
        return dataframe

    def _request_boqa(self, hpo_ids: [str], prefilter: {str: str}={}) \
            -> pandas.DataFrame:
        '''Get boqa information on a list of hpo ids.
        '''
        hpo_ids = ','.join(hpo_ids)
        boqa_names = ['value', 'nothing', 'disease-id', 'disease-name']
        params = {
            'username': self.user,
            'password': self.password,
            'mobilequery': 'true',
            'doboqa': 'true',
            'terms': hpo_ids,
            'numres': 100
                }
        dataframe = self._request_data_as_df(
            params=params, names=boqa_names, prefilter=prefilter)
        # remove unneeded columns, the nothing column is literally nothing
        dataframe = dataframe.drop(['nothing'], axis=1)
        return dataframe
