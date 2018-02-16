'''
OMIM Api Bindings
---
This class handles automatic download of morbidmap and mim2gene.

In the future specific usage of the OMIM API might be implemented.
'''
import re
import os
import logging

import pandas
import requests

RE_OMIM_PHEN = re.compile('.* (\d{6}) \((\d)\)')


class Omim:
    '''
    Download omim files from official website resources or get them locally
    if they already exist
    '''

    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) \
        AppleWebKit/537.36 (KHTML, like Gecko) \
        Chrome/39.0.2171.95 Safari/537.36'
    }

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
        if config:
            self._use_cached = config.omim['use_cached']
            self._api_key = config.omim['api_key']
            self._mimdir = config.omim['mimdir']
        else:
            self._use_cached = use_cached
            self._api_key = api_key
            self._mimdir = mimdir

        # create directory if it doesnt exist
        os.makedirs(mimdir, exist_ok=True)

        mim2gene = self._retrieve_file(
            os.path.join(mimdir, 'mim2gene.txt'),
            'https://omim.org/static/omim/data/mim2gene.txt',
            ['mim_number', 'mim_entry_type', 'entrez_id', 'gene_symbol',
             'ensembl'])
        self.mim2gene = mim2gene.dropna(
            subset=['mim_number', 'entrez_id']).set_index('mim_number')
        self.entrez2gene = mim2gene.dropna(
            subset=['mim_number', 'entrez_id']).set_index('entrez_id')

        morbidmap_url = \
            'https://data.omim.org/downloads/{}/morbidmap.txt'.format(api_key)
        morbidmap = self._retrieve_file(
            os.path.join(mimdir, 'morbidmap.txt'), morbidmap_url,
            ['phenotype', 'gene_symbol', 'mim_number', 'cyto_location'])
        morbidmap['phen_mim_number'] = \
            morbidmap['phenotype'].apply(self._extract_omim)
        self.morbidmap = morbidmap

        self.phen_to_mim = self._create_phen_to_mim()

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
        if not os.path.exists(path) or not self._use_cached:
            resp = requests.get(url, headers=self.headers)
            # save the downloaded file to a disk file
            with open(path, 'wb') as savedfile:
                for data in resp.iter_content():
                    savedfile.write(data)
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
            logging.warning("%s is ambiguous", str(needle))
        return entry[target]

    def mim_gene_to_entrez_id(self, mim_gene):
        try:
            entrez_id = self.mim2gene.at[str(mim_gene), 'entrez_id']
        except KeyError:
            entrez_id = ''
        return entrez_id

    def mim_gene_to_symbol(self, mim_gene):
        gene_symbol = self.mim2gene.at[str(mim_gene), 'gene_symbol']
        return gene_symbol

    def entrez_id_to_mim_gene(self, entrez_id):
        try:
            gene_omim_id = self.entrez2gene.at[str(entrez_id), 'mim_number']
        except KeyError:
            gene_omim_id = ''
        return gene_omim_id

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


def main():
    '''This should only be used for testing purposes.
    '''
    omim = Omim(api_key='<INSERT CORRECT API KEY>', mimdir='test',
                use_cached=False)
    print(omim.morbidmap)


if __name__ == '__main__':
    main()
