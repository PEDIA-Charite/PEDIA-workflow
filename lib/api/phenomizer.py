import io
import logging
import re

import requests
import pandas

'''
Phenomization
---
Using a list of HPO features a list of associated diseases with causal genes and a likelyhood score
is returned.

Because our current approach is on a gene-level scoring. This information will be rearranged into a
gene-to-score mapping. Thus the associated diseases returned by phenomization and boqa can be different.

'''

RE_SYMBOL = re.compile('(\w+) \(\d+\)')
def match_symbol(string):
    m = RE_SYMBOL.search(string)
    return m and m.group(1) or string

def get_max_gene(df_group):
    i = df_group['value'].idxmax()
    row = df_group.loc[i]
    return row

class PhenomizerService(requests.Session):
    '''Handling of interop with Phenomizer service, which provides the pheno and boqa scores used in the process.
    '''

    def __init__(self, Phenomizer_Url='', Phenomizer_User='', Phenomizer_Password='', config=None):
        '''
        Create a new phenomizer service instance.

        Params:
            url: Url of phenomizer service
        '''
        super().__init__()
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
                total = 3
                ,read = 3
                ,connect = 3
                ,backoff_factor = 0.3
                , status_forcelist=(500,)
                )
        adapter = requests.adapters.HTTPAdapter(max_retries = retry)

        self.mount('http://', adapter)
        self.mount('https://', adapter)

    def get_df(self, url, params, names, prefilter):
        '''prefilter - Set whether dataframes should be filtered down to disease-id containing OMIM.
        '''
        r = self.get(url, params=params)
        r.raise_for_status()
        rstring = "\n".join([ s for s in r.text.split('\n') if not s.startswith('#') or s == '' ])
        rawdata = io.StringIO(rstring)
        df = pandas.read_table(rawdata, sep='\t', index_col=None, header=None, names=names)
        if prefilter:
            for k,v in prefilter.items():
                df = df.loc[df[k].str.contains(v)]
        return df

    def disease_boqa_phenomize(self, hpo_ids):
        if not hpo_ids:
            return None
        # this process might need to be retried, but it is currently reliable enough to run directly
        hpo_df = self._request_phenomize(hpo_ids, prefilter={'disease-id':'OMIM'})
        boqa_df = self._request_boqa(hpo_ids, prefilter={'disease-id':'OMIM'})

        # preprocess result for join based on omim disease id
        hpo_df['value'] = 1 - hpo_df['value']
        # remove id tags
        hpo_df['disease-id'] = hpo_df['disease-id'].apply(lambda x: x.split(':')[-1])
        hpo_df = hpo_df.set_index('disease-id')
        boqa_df['disease-id'] = boqa_df['disease-id'].apply(lambda x: x.split(':')[-1])
        boqa_df = boqa_df.set_index('disease-id')
        # join dataframes on disease id
        scores_df = hpo_df.join(boqa_df, how='outer', lsuffix='_pheno', rsuffix='_boqa')
        scores_df = scores_df.fillna(
                {'value_pheno':0.0, 'value_boqa':0.0,
                    'gene-symbol':'', 'gene-id': '',
                    'disease-name_boqa' : '', 'disease-name_pheno': ''}
                )
        return scores_df

    def _request_phenomize(self, hpo_ids, prefilter={}):
        hpo_ids = ','.join(hpo_ids)
        phen_names = ['value','score','disease-id','disease-name','gene-symbol','gene-id']
        params = {
                'mobilequery' : 'true'
                ,'username' : self.user
                ,'password' : self.password
                ,'terms' : hpo_ids
                ,'numres' : 100
                }
        df = self.get_df(self.url, params=params,names=phen_names, prefilter=prefilter)
        df = df.drop(['score'], axis=1)
        return df

    def _request_boqa(self, hpo_ids, prefilter={}):
        hpo_ids = ','.join(hpo_ids)
        boqa_names = ['value', 'nothing', 'disease-id','disease-name']
        params = {
                'username' : self.user
                ,'password' : self.password
                ,'mobilequery' : 'true'
                ,'doboqa' : 'true'
                ,'terms' : hpo_ids
                ,'numres' : 100
                }
        df = self.get_df(self.url, params=params,names=boqa_names, prefilter=prefilter)
        df = df.drop(['nothing'], axis=1)
        return df
