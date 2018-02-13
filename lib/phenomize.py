import json
import requests
import os
import sys
import re
import io
from argparse import ArgumentParser
import logging
import configparser

import pandas
import numpy

'''
jSON Phenomization tools
---
These functions add phenomization scores from Phenomizer and Boqa to the case information available.
Phenomization and Boqa entries without genetic correlate will also be added back to the list of genetic candidates.
'''

class Annotator:
    '''Add phenomization and boqa scores to the geneList present in the json file.
    '''

    def __init__(self, pheno_args, omim_args):
        self.phenomizer = PhenomizerService(**pheno_args)
        self.omim = get_omim_files(api_key=omim_args['key'])


    def fill_missing(self, row):
        omim = self.get_omim_with_id(row['gene_id'])
        if omim.empty:
            return row
        omim = omim.iloc[0]
        if omim.empty:
            return row
        if pandas.isna(row['gene_omim_id']):
            row['gene_omim_id'] = omim['mim_number']
        if pandas.isna(row['gene_symbol']):
            row['gene_symbol'] = omim['gene_symbol']
        logging.info(
                "Filled missing for gene_id {} with {} mim and {} symbol".format(
                    row['gene_id'],row['gene_omim_id'],row['gene_symbol'])
                )
        return row

    def process(self, inputjs):
        outputjs = {}
        outputjs = inputjs
        outputjs['processing']=['python phenomize_jsons.py']

        bq = self.get_boqa(inputjs)
        bq_df = pandas.DataFrame({'gene_id':list(bq.keys()),'boqa_score':list(bq.values())})
        pheno = self.get_pheno(inputjs)
        pheno_df = pandas.DataFrame({'gene_id':list(pheno.keys()),'pheno_score':list(pheno.values())})

        scores = bq_df.merge(pheno_df, how='outer')
        entries = pandas.DataFrame(inputjs['geneList'])
        new_data = entries.merge(scores, left_on='gene_id', right_on='gene_id',how='outer')
        # we still need omim id and gene symbol for all entries
        new_data = new_data.apply(self.fill_missing,axis=1)

        new_dicts = new_data.to_dict('records')

        outputjs['geneList'] = new_dicts
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
    '''Generator class to serve json files from the specified folder and handle saving of edited objects back as jsons.
    '''

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

    def save(self, filename, obj):
        filepath = os.path.join(self.outputdir, filename)
        with open(filepath, 'w') as f:
            json.dump(obj, f)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('inputdir')
    parser.add_argument('outputdir')
    parser.add_argument('--configfile','-c', type=str, default='config.ini')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.configfile)
    phen_config = config['phenomization']

    omim_args = {'key':phen_config['OMIM_Key'],'cached':phen_config['OMIM_Cached']}
    phenom_args = {'url':phen_config['Phenomizer_Url']
            ,'user':phen_config['Phenomizer_User']
            ,'password':phen_config['Phenomizer_Password']}

    annot = Annotator(pheno_args=phenom_args,omim_args=omim_args)

    files = JsonProcessor(args.inputdir, args.outputdir)

    for filename, js in files:
        annot.process(js)
        files.save(filename, annot.process(js))
