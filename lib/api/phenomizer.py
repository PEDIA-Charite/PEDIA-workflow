import io
import logging
import re

import requests
import pandas

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

    phen_names = ['value','score','disease-id','disease-name','gene-symbol','gene-id']
    boqa_names = ['value', 'nothing', 'disease-id','disease-name']
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

    def get_df(self, url, params, names):
        r = self.get(url, params=params)
        r.raise_for_status()
        rstring = "\n".join([ s for s in r.text.split('\n') if not s.startswith('#') or s == '' ])
        rawdata = io.StringIO(rstring)
        df = pandas.read_table(rawdata, sep='\t', index_col=None, header=None, names=names)
        return df

    def pheno_boqa(self, hpo_ids, omim):
        pheno = self.get_phenomize(hpo_ids)
        boqa = self.get_boqa(hpo_ids, omim)
        scores = pheno.join(boqa, on='gene-id', how='outer', lsuffix='_pheno', rsuffix='_boqa')
        scores['gene-symbol'] = scores[['gene-symbol_pheno','gene-symbol_boqa']].apply(lambda x: pandas.isna(x[1]) and x[0] or x[1], axis=1)
        scores = scores.drop(['gene-symbol_boqa', 'gene-symbol_pheno', 'gene-id_boqa', 'gene-id_pheno'], axis=1)
        scores = scores.fillna({'value_boqa' : 0.0, 'value_pheno' : 0.0})
        scores = scores.reset_index(drop=True)
        return scores

    def _request_phenomize(self, hpo_ids):
        params = {
                'mobilequery' : 'true'
                ,'username' : self.user
                ,'password' : self.password
                ,'terms' : hpo_ids
                ,'numres' : 100
                }
        df = self.get_df(self.url, params=params,names=self.phen_names)
        return df

    def _request_boqa(self, hpo_ids):
        params = {
                'username' : self.user
                ,'password' : self.password
                ,'mobilequery' : 'true'
                ,'doboqa' : 'true'
                ,'terms' : hpo_ids
                ,'numres' : 100
                }
        df = self.get_df(self.url, params=params,names=self.boqa_names)
        return df

    def get_phenomize(self, hpo_ids):
        if hpo_ids == '':
            return {}
        try:
            hpo_df = self._request_phenomize(hpo_ids)
        except:
            logging.warning("{} generated request error".format(hpo_ids))
            return {}
        hpo_df = hpo_df.drop(['score', 'disease-name'], axis=1)
        hpo_df['value'] = 1 - hpo_df['value']
        hpo_df = hpo_df.dropna(subset=['gene-id'])
        s = hpo_df.apply(lambda x: pandas.Series(x['gene-id'].split(', ')), axis=1).stack().reset_index(level=1, drop=True)
        symbols = hpo_df.apply(lambda x: pandas.Series(x['gene-symbol'].split(', ')), axis=1).stack().reset_index(level=1, drop=True)
        s.name = 'gene-id'
        hpo_df = hpo_df.drop(['gene-id'],axis=1).join(s)
        hpo_df = hpo_df.drop(['gene-symbol'],axis=1)
        hpo_df['gene-symbol'] = symbols
        results = hpo_df.groupby(['gene-id']).apply(get_max_gene)
        results['gene-symbol'] = results['gene-symbol'].apply(match_symbol)
        results['disease-id'] = results['disease-id'].apply(lambda x: x.split(':',1)[1])
        return results

    def get_boqa(self, hpo_ids, omim):
        if hpo_ids == '':
            return {}
        try:
            boqa_df = self._request_boqa(hpo_ids)
        except:
            logging.warning("{} generated boqa request error".format(hpo_ids))
        boqa_df = boqa_df.drop(['nothing', 'disease-name'], axis=1)
        boqa_df = boqa_df.loc[boqa_df['disease-id'].str.contains('OMIM')]
        boqa_df['disease-id'] = boqa_df['disease-id'].apply(lambda x: x.split(':')[2])

        t = boqa_df['disease-id'].apply(lambda x: omim.mim_pheno_to_gene(x))
        t = t.loc[t.apply(lambda x: bool(x))]
        boqa_df['res'] = t
        boqa_df = boqa_df.dropna(subset=['res'])
        s = boqa_df.apply(lambda x: pandas.Series(x['res']), axis=1).stack().reset_index(level=1, drop=True)
        s.name = 'res'
        stacked = boqa_df.drop(['res'],axis=1).join(s)

        stacked['gene-id'] = stacked['res'].apply(lambda x: x['gene_id'])
        stacked['gene-symbol'] = stacked['res'].apply(lambda x: x['gene_symbol'])
        results = stacked.groupby(['gene-id']).apply(get_max_gene)
        results = results.drop(['res'], axis=1)
        return results
