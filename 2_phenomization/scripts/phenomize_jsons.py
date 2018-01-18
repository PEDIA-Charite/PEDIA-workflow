import json
import requests
import os
import sys
import re
from argparse import ArgumentParser
import requests
import logging

import pandas
import numpy
import io

'''
jSON Phenomization tools
---

'''

RE_OMIM_PHEN = re.compile('.* (\d{6}) \((\d)\)')

def retrieve_file(path, url, cached, names):
    '''
    Get file from local filesystem. If it doesnt exists or cached is false, save from remote url to the local path
    and thereafter read from local filesystem.
    '''
    if not os.path.exists(path) or not cached:
        r = requests.get(url)
        with open(path, 'wb') as f:
            for data in r.iter_content():
                f.write(data)
    return pandas.read_table(path, delimiter='\t', comment='#', names=names, dtype=str)


def extract_omim(raw_omim):
    if pandas.isna(raw_omim):
        return numpy.nan
    match = RE_OMIM_PHEN.search(raw_omim)
    if match is None:
        return numpy.nan
    return match.group(1)

def get_omim_files(mimdir='', api_key='', use_cached=True):
    '''
    Download omim files from official website resources or get them locally if they already exist

    Args:
        mimdir: Directory where omim files are saved and retrieved, defaults to the current working directory
        api_key: OMIM API key retrieve it from the official omim website
        use_cached: Whether already downloaded files are used, if False omim files are always downloaded

    Returns:
        Data stucture with mim2gene and morbidmap
    '''
    mim2gene = retrieve_file('mim2gene.txt'
            ,'https://omim.org/static/omim/data/mim2gene.txt'
            ,use_cached
            ,['mim_number','mim_entry_type','entrez_id', 'gene_symbol','ensembl'])

    morbidmap_url = 'https://data.omim.org/downloads/{}/morbidmap.txt'.format(api_key)
    morbidmap = retrieve_file('morbidmap.txt', morbidmap_url, use_cached
            ,['phenotype', 'gene_symbol', 'mim_number', 'cyto_location'])

    morbidmap['phen_mim_number'] = morbidmap['phenotype'].apply(extract_omim)

    return {
            'mim2gene' : mim2gene
            ,'morbidmap' : morbidmap
            }

class PhenomizerService(requests.Session):
    phen_names = ['value','score','disease-id','disease-name','gene-symbol','gene-id']
    boqa_names = ['p', 'nothing', 'disease-id','disease-name']
    def __init__(self, url, user, password):
        '''
        Create a new phenomizer service instance.

        Params:
            url: Url of phenomizer service
        '''
        super().__init__()
        self.url = url
        self.user = user
        self.password = password
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

    def request_phenomize(self, hpo_ids):
        params = {
                'mobilequery' : 'true'
                ,'username' : self.user
                ,'password' : self.password
                ,'terms' : hpo_ids
                ,'numres' : 100
                }
        df = self.get_df(self.url, params=params,names=self.phen_names)
        return df

    def request_boqa(self, hpo_ids):
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

class Annotator:
    def __init__(self, pheno_args, omim_args):
        self.phenomizer = PhenomizerService(**pheno_args)
        self.omim = get_omim_files(api_key=omim_args['key'])


    def get_morbidmap(self, omim_id):
        mm_entry = self.omim['morbidmap'].loc[self.omim['morbidmap']['phen_mim_number'] == omim_id]
        return mm_entry

    def get_omim_with_id(self, gene_id):
        omim_entry = self.omim['mim2gene'].loc[self.omim['mim2gene']['entrez_id'].astype('float32') == float(gene_id)]
        return omim_entry

    def get_omim_with_name(self, gene_name):
        omim_entry = self.omim['mim2gene'].loc[self.omim['mim2gene']['gene_symbol'] == gene_name]
        return omim_entry

    def process(self, inputjs):
        outputjs = {}
        outputjs = inputjs
        outputjs['processing']=['python phenomize_jsons.py']

        bq = self.get_boqa(inputjs)
        bq_df = pandas.DataFrame({'score_id':list(bq.keys()),'score_boqa':list(bq.values())})
        pheno = self.get_pheno(inputjs)
        pheno_df = pandas.DataFrame({'score_id':list(pheno.keys()),'score_pheno':list(pheno.values())})

        scores = bq_df.merge(pheno_df, how='outer')
        entries = pandas.DataFrame(inputjs['geneList'])
        new_data = entries.merge(scores, left_on='gene_id', right_on='score_id',how='outer')

        outputjs['geneList'] = new_data.to_dict('records').values()
        return outputjs

    def get_pheno(self, inputjs):
        hpo_ids = ",".join(inputjs['features'])
        if hpo_ids == '':
            return {}
        try:
            hpo_df = self.phenomizer.request_phenomize(hpo_ids)
        except:
            logging.warning("{} generated request error".format(hpo_ids))
            return {}
        hpo_df['value'] = 1 - hpo_df['value']
        score_dict = {}
        for i,h in hpo_df.iterrows():
            if pandas.isna(h['gene-symbol']):
                continue
            gene_ids = [ x.strip() for x in h['gene-id'].split(',') ]
            genes = [ x.strip() for x in h['gene-symbol'].split(',') ]
            assert len(gene_ids) == len(genes) ,'{} gene ids and {} genes'.format(gene_ids,genes)
            ## unnecessary as long as we assume that all genes have ids
            # genes = [ RE_GENE_SYMBOL.match(x)
            #         for x in genes if RE_GENE_SYMBOL.match(x) is not None ]
            for g in gene_ids:
                if g not in score_dict:
                    score_dict[g] = h['value']
                else:
                    score_dict[g] = max(h['value'],score_dict[g])

        return score_dict

    def get_boqa(self, inputjs):
        hpo_ids = ",".join(inputjs['features'])
        if hpo_ids == '':
            return {}
        try:
            boqa_df = self.phenomizer.request_boqa(hpo_ids)
        except:
            logging.warning("{} generated boqa request error".format(hpo_ids))
        score_dict = {}
        for i,h in boqa_df.iterrows():
            if pandas.isna(h['disease-id']):
                continue
            if 'OMIM' not in h['disease-id']:
                continue
            mim = h['disease-id'].split(':')[2]
            mm_genes = self.get_morbidmap(mim)
            if mm_genes.empty:
                logging.warning('OMIM {} not contained in morbidmap.'.format(mim))
                continue
            if mm_genes.shape[0] > 1:
                mm_genes = [ y.split(',') for y in mm_genes['gene_symbol'] ]
                mm_genes = [ x.strip() for y in mm_genes for x in y ]
            else:
                mm_genes = mm_genes['gene_symbol']
                mm_genes = mm_genes.iloc[0]
                mm_genes = [ x.strip() for x in mm_genes.split(',') ]
            gene_ids = [ self.get_omim_with_name(x)
                for x in mm_genes ]
            gene_ids = [ x['entrez_id'].iloc[0] for x in gene_ids if not x.empty ]
            for gi in gene_ids:
                if gi not in score_dict:
                    score_dict[gi] = h['p']
                else:
                    score_dict[gi] = max(h['p'],score_dict[gi])
        return score_dict



class JsonProcessor:
    def __init__(self, inputdir, outputdir):
        self._cache = {}
        self.inputdir = inputdir.rstrip('/\\')
        self.inputindex = 0
        self.outputdir = outputdir.rstrip('/\\')

        if not os.path.exists(self.outputdir):
            os.makedirs(self.outputdir)

        self.filenames = os.listdir(self.inputdir)

    def __iter__(self):
        return self

    def __next__(self):
        if self.inputindex < len(self.filenames) - 1:
            cur, self.inputindex = self.inputindex, self.inputindex + 1
            return self.load(self.filenames[cur])
        else:
            raise StopIteration()

    def load(self, fname):
        fobj = json.load(open(os.path.join(self.inputdir,fname), 'r'))
        return fname, fobj

    def save(self, filename, json):
        filepath = os.path.join(self.outputdir, filename)
        with open(filepath, 'w') as f:
            json.dump(json, f)

def phenomize_argv():
    parser = ArgumentParser()
    parser.add_argument('inputdir')
    parser.add_argument('outputdir')
    args = parser.parse_args()
    return args.inputdir, args.outputdir


if __name__ == '__main__':
    phenomizer_args = {
            'url' : os.getenv('PHENOMIZER_URL')
            ,'user' : os.getenv('PHENOMIZER_USER')
            ,'password' : os.getenv('PHENOMIZER_PW')
    }
    omim_args = {
            'key' : os.getenv('OMIM_KEY')
            ,'cached' : os.getenv('OMIM_CACHED').upper() == 'TRUE'
    }

    if omim_args['key'] is None and not omim_args['cached']:
        print(
'''ease set OMIM_KEY in your environment variables to retrieve OMIM morbidmap.
Altnatively provide genemap and morbidmap in the script directory and set OMIM_CACHED to TRUE'''
        )
        sys.exit(1)

    annot = Annotator(pheno_args=phenomizer_args,omim_args=omim_args)

    files = JsonProcessor(*phenomize_argv())

    for filename, js in files:
        annot.process(js)
        files.save(filename, annot.process(js))
