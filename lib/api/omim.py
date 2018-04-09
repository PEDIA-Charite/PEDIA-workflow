'''
OMIM Api Bindings
---
This class handles automatic download of morbidmap and mim2gene.

In the future specific usage of the OMIM API might be implemented.
'''
import re
import os
import json
import logging

import pandas
import requests

from lib.utils import get_file_hash

RE_OMIM_PHEN = re.compile('.* (\d{6}) \((\d)\)')

LOGGER = logging.getLogger(__name__)


class Omim:
    '''
    Download omim files from official website resources or get them locally
    if they already exist
    '''

    url = "https://api.omim.org/api/"

    def __init__(self, api_key: str='', mimdir: str='data',
                 use_cached: bool=True, config: 'ConfigManager'=None):
        '''
        Omim API requires an API key for the download of the morbidmap.

        Args:
            mimdir: Directory where omim files are saved and retrieved,
                    defaults to the current working directory

            api_key: OMIM API key retrieve it from the official omim website

            use_cached: Whether already downloaded files are used, if
                        False omim files are always downloaded

        Returns:
            Data stucture with mim2gene and morbidmap
        '''
        # use config object to get parameters
        self._use_cached = use_cached
        self._api_key = api_key
        self._mimdir = mimdir

        mim2gene_hash = ''
        morbidmap_hash = ''

        if config:
            self._use_cached = config.omim['use_cached']
            self._api_key = config.omim['api_key']
            self._mimdir = config.omim['mimdir']
            mim2gene_hash = config.omim['mim2gene_hash']
            morbidmap_hash = config.omim['morbidmap_hash']

        # save hashes to object for usage as stamp
        self.mim2gene_has = mim2gene_hash
        self.morbidmap_hash = morbidmap_hash

        # create directory if it doesnt exist
        os.makedirs(mimdir, exist_ok=True)

        # add api key to http header
        self.headers = {
            "ApiKey": self._api_key
        }

        mim2gene_path = os.path.join(mimdir, 'mim2gene.txt')
        # ensure using revision specified in config
        if mim2gene_hash:
            assert get_file_hash(mim2gene_path) == mim2gene_hash

        mim2gene = self._retrieve_file(
            mim2gene_path,
            'https://omim.org/static/omim/data/mim2gene.txt',
            ['mim_number', 'mim_entry_type', 'entrez_id', 'gene_symbol',
             'ensembl'])

        self.mim2gene = mim2gene.dropna(
            subset=['mim_number', 'entrez_id']
        ).drop_duplicates(subset="entrez_id").set_index('mim_number')

        self.entrez2gene = mim2gene.dropna(
            subset=['mim_number', 'entrez_id']
        ).drop_duplicates(subset="entrez_id").set_index('entrez_id')

        morbidmap_url = \
            'https://data.omim.org/downloads/{}/morbidmap.txt'.format(api_key)

        morbidmap_path = os.path.join(mimdir, 'morbidmap.txt')
        # ensure using revision specified in config
        if morbidmap_hash:
            assert get_file_hash(morbidmap_path) == morbidmap_hash

        morbidmap = self._retrieve_file(
            morbidmap_path, morbidmap_url,
            ['phenotype', 'gene_symbol', 'mim_number', 'cyto_location'])
        morbidmap['phen_mim_number'] = \
            morbidmap['phenotype'].apply(self._extract_omim)
        self.morbidmap = morbidmap

        self.phen_to_mim = self._create_phen_to_mim()

        # load mim_to_ps json file
        self.mim_to_ps = {}
        filepath = os.path.join(self._mimdir, "mim_to_ps.json")
        if os.path.exists(filepath):
            with open(filepath, "r") as jsfile:
                jsdata = json.load(jsfile)
            # check checksum of the mapping file
            if jsdata["morbidmap_hash"] == self.morbidmap_hash:
                self.mim_to_ps = jsdata["mapping"]
            else:
                LOGGER.warning(
                    ("mim_to_ps.json MD5 hash for morbidmap does "
                     "not match current morbidmap checksum.")
                )

    def _create_phen_to_mim(self):
        '''Create dictionary mapping phenotypic omim ids to list of gene ids.
        '''
        pheno_entries = self.morbidmap.to_dict('records')
        phen_to_mim = {}
        for entry in pheno_entries:
            phen_mim = entry['phen_mim_number']
            if phen_mim is None:
                continue
            phen_mim = str(phen_mim)
            entry = (entry['mim_number'], entry['gene_symbol'].split(', ')[0])
            if phen_mim in phen_to_mim:
                phen_to_mim[phen_mim].append(entry)
            else:
                phen_to_mim[phen_mim] = [entry]
        return phen_to_mim

    def _retrieve_file(self, path, url, names):
        '''Get file from local filesystem.
        If it doesnt exists or cached is false, save from remote url to the
        local path and thereafter read from local filesystem.
        '''
        # if not os.path.exists(path) or not self._use_cached:
        #     resp = requests.get(url, headers=self.headers)
        #     # save the downloaded file to a disk file
        #     with open(path, 'wb') as savedfile:
        #         for data in resp.iter_content():
        #             savedfile.write(data)
        if not os.path.exists(path):
            raise RuntimeError(
                "{} does not exist. Add them to the location manually.")
        return pandas.read_table(
            path, delimiter='\t', comment='#', names=names, dtype=str)

    def _extract_omim(self, raw_omim):
        '''Extract phenotypic MIM number from the phenotypic string.
        '''
        if pandas.isna(raw_omim):
            return None
        match = RE_OMIM_PHEN.search(raw_omim)
        if match is None:
            return None
        return match.group(1)

    def _search_single(self, dataframe, column, needle, target):
        '''Search a single entry in a dataframe on a per-row basis.
        '''
        # pandas saves values as object type per default, coerce both to str to
        # make comparisons easier to debug correctly
        entry = dataframe.loc[dataframe[column] == str(needle)]
        if entry.empty:
            return []
        # unpack the dataframe to get the single cell of interest
        # tell us if our requests are returning multiple entries
        if entry.shape[0] > 1:
            LOGGER.warning("%s is ambiguous", str(needle))
        return entry[target]

    def search_table(self, table, needle, target, default=""):
        try:
            result = table.at[needle, target]
        except KeyError:
            result = default
        return result

    def mim_gene_to_entrez_id(self, mim_gene):
        return self.search_table(
            self.mim2gene,
            str(mim_gene),
            "entrez_id"
        )

    def mim_gene_to_symbol(self, mim_gene):
        return self.search_table(
            self.mim2gene,
            str(mim_gene),
            "gene_symbol"
        )

    def entrez_id_to_mim_gene(self, entrez_id):
        return self.search_table(
            self.entrez2gene,
            str(entrez_id),
            "mim_number"
        )

    def entrez_id_to_symbol(self, entrez_id):
        return self.search_table(
            self.entrez2gene,
            str(entrez_id),
            "gene_symbol"
        )

    def mim_pheno_to_mim_gene(self, mim_pheno):
        gene_omim_id = self._search_single(
            self.morbidmap, 'phen_mim_number', mim_pheno, 'mim_number')
        return gene_omim_id

    def mim_pheno_to_gene(self, mim_pheno):
        '''Convert a phenotype omim id to a gene dictionary object containing
        keys: gene_id, gene_symbol and gene_omim_id.
        '''
        if str(mim_pheno) in self.phen_to_mim:
            gene_entry = {}
            for mim_gene, gene_symbol in self.phen_to_mim[str(mim_pheno)]:
                entrez_id = self.mim_gene_to_entrez_id(mim_gene)
                gene_entry[mim_gene] = {'gene_id': entrez_id,
                                        'gene_symbol': gene_symbol,
                                        'gene_omim_id': mim_gene}
            return gene_entry
        else:
            return {}

    def api_query_entry(self, mim_pheno: str):
        '''Get entry information for a phenotypic omim id.'''
        params = {
            'format': 'json',
            'include': 'existFlags,geneMap',
            'mimNumber': mim_pheno
        }
        url = self.url + "entry"
        resp = requests.get(url, params=params, headers=self.headers)
        js_data = [
            e["entry"]
            for e in resp.json()["omim"]["entryList"]
        ]
        return js_data

    def api_query_phenotype_mapping(self, mim_pheno: str):
        '''Map phenotype to phenotypic series number.'''
        entry_list = self.api_query_entry(mim_pheno)
        phenotypic_numbers = [
            v['phenotypeMap']['phenotypicSeriesNumber']
            for e in entry_list
            if 'phenotypeMapList' in e
            for v in e['phenotypeMapList']
            if 'phenotypicSeriesExists' in e and e['phenotypicSeriesExists']
        ]
        return phenotypic_numbers

    def construct_phenotypic_series_mapping(self):
        '''Construct a dict with mapping of phenotypic omim ids to
        phenotypic series ids.'''

        output_path = os.path.join(self._mimdir, "mim_to_ps.json")
        if os.path.exists(output_path):
            with open(output_path, "r") as old_json:
                old = json.load(old_json)
                if old['morbidmap_hash'] == self.morbidmap_hash:
                    LOGGER.warning(
                        ("Phenotypic series mapping with same base "
                         "morbidmap has already been constructed.")
                    )
                    return
        mim_numbers = self.morbidmap[
            "phen_mim_number"].dropna().drop_duplicates()

        mim_ps_mapping = {
            num: self.api_query_phenotype_mapping(num)
            for num in mim_numbers
        }
        mim_ps_data = {
            'morbidmap_hash': self.morbidmap_hash,
            'mapping': mim_ps_mapping
        }
        with open(output_path, "w") as new_json:
            json.dump(mim_ps_data, new_json)
        return mim_ps_data

    def omim_id_to_phenotypic_series(self, omim_id: str) -> str:
        '''Translate omim id to phenotypic series id, or if it does
        not exist to the empty string.'''
        assert isinstance(omim_id, str), "OMIM ID has to be string."

        ps_label = ""
        # take the first phenotypic series entry or if empty leave the
        # empty string
        if omim_id in self.mim_to_ps:
            ps_label = self.mim_to_ps[omim_id][0] \
                if self.mim_to_ps[omim_id] else ps_label
        return ps_label


def main():
    '''This should only be used for testing purposes.
    '''
    omim = Omim(api_key='<INSERT CORRECT API KEY>', mimdir='test',
                use_cached=False)
    print(omim.morbidmap)


if __name__ == '__main__':
    main()
