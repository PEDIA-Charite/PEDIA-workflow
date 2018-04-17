'''
OMIM Api Bindings
---
This class handles automatic download of morbidmap and mim2gene.

In the future specific usage of the OMIM API might be implemented.
'''
import re
import os
import csv
import json
import logging
import collections
import typing

import pandas

from lib.utils import get_file_hash

RE_OMIM_PHEN = re.compile(r'.* (\d{6}) \((\d)\)')

LOGGER = logging.getLogger(__name__)


def id_to_string(val):
    if isinstance(val, int):
        val = str(val)
    elif isinstance(val, float):
        LOGGER.warning(
            "Conversion of val from float value %f", val
        )
        val = str(int(val))
    return val


def omim_check(func):
    '''Coerce omim ids to string type.'''

    def checker(obj, omim_id, *args, **kwargs):
        omim_id = id_to_string(omim_id)
        if not isinstance(omim_id, str):
            raise TypeError("ID not string", omim_id)
        return func(obj, omim_id, *args, **kwargs)

    return checker


def omim_list(func):
    '''Coerce list or single data to list format.'''

    def checker(obj, omim_ids, *args, **kwargs):
        omim_ids = id_to_string(omim_ids)
        if isinstance(omim_ids, str):
            omim_ids = [omim_ids]
        elif omim_ids is None:
            omim_ids = []
        return func(obj, omim_ids, *args, **kwargs)
    return checker


class Omim:
    '''
    Download omim files from official website resources or get them locally
    if they already exist
    '''

    url = "https://api.omim.org/api/"

    data_meta = {
        "morbidmap": {
            "filename": "morbidmap.txt",
            "type": "pandas",
            "colnames": [
                'phenotype', 'gene_symbol', 'mim_number', 'cyto_location'
            ]
        },
        "mim2gene": {
            "filename": "mim2gene.txt",
            "type": "pandas",
            "colnames": [
                'mim_number', 'mim_entry_type',
                'entrez_id', 'gene_symbol', 'ensembl'
            ]
        },
        "phenotypicSeries": {
            "filename": "phenotypicSeries.txt",
            "type": "ps_data",
            "colnames": [
                "ps_number", "mim_number", "phenotype"
            ]
        },
        "omim_deprecated_replacement": {
            "filename": "omim_deprecated_replacement.json",
            "type": "json",
            "colnames": []
        },
    }

    indexes_meta = {
        "gene_omim": {
            "source": "mim2gene"
        },
        "entrez_id": {
            "source": "mim2gene"
        },
        "pheno_omim": {
            "source": "morbidmap"
        },
        "pheno_series": {
            "source": "phenotypicSeries"
        }
    }

    def __init__(
            self,
            mimdir: str = 'data',
            config: 'ConfigManager' = None
    ):
        '''
        Omim API requires an API key for the download of the morbidmap.

        Args:
            mimdir: Directory where omim files are saved and retrieved,
                    defaults to the current working directory

            api_key: OMIM API key retrieve it from the official omim website

        Returns:
            Data stucture with mim2gene and morbidmap
        '''
        # use config object to get parameters
        self._mimdir = mimdir

        mim2gene_hash = ''
        morbidmap_hash = ''

        if config:
            self._mimdir = config.omim['mimdir']
            mim2gene_hash = config.omim['mim2gene_hash']
            morbidmap_hash = config.omim['morbidmap_hash']

        # save hashes to object for usage as stamp
        hashes = {
            "mim2gene": mim2gene_hash,
            "morbidmap": morbidmap_hash
        }

        # create directory if it doesnt exist
        os.makedirs(mimdir, exist_ok=True)

        self.files = {}
        for filename, fileinfo in self.data_meta.items():
            filehash = hashes[filename] if filename in hashes else None
            self.files[filename] = self.load_file(mimdir, fileinfo, filehash)
            self.files[filename] = self.post_ops(
                self.files[filename], filename
            )

        self.indexes = {}
        for index_name, index_info in self.indexes_meta.items():
            self.indexes[index_name] = self.create_index(
                index_name, index_info
            )

    def create_index(self, name, info):
        '''Create index for faster access.'''
        data = self.files[info["source"]]
        if name == "entrez_id":
            index = data.dropna(
                subset=['mim_number', 'entrez_id']
            ).drop_duplicates(subset="entrez_id").set_index('entrez_id')
        elif name == "gene_omim":
            index = data.dropna(
                subset=['mim_number', 'entrez_id']
            ).drop_duplicates(subset="entrez_id").set_index('mim_number')
        elif name == "pheno_omim":
            index = self.create_phen_to_mim(data)
        elif name == "pheno_series":
            index = self.create_phen_omim_to_ps(data)
        return index

    def load_file(self, mimdir, fileinfo, filehash):
        '''Load file depending on type with optional hash checking.'''
        filepath = os.path.join(mimdir, fileinfo["filename"])
        if filehash:
            assert get_file_hash(filepath) == filehash

        if fileinfo["type"] == "ps_data":
            data = self.load_phenotypic_series(
                filepath, fileinfo["colnames"]
            )
        elif fileinfo["type"] == "pandas":
            data = self.load_dataframe(filepath, fileinfo["colnames"])
        elif fileinfo["type"] == "json":
            with open(filepath, "r") as handle:
                data = json.load(handle)
        return data

    def post_ops(self, data, filename):
        '''Post processing operations depending on data name.'''
        if filename == "morbidmap":
            data['phen_mim_number'] = data['phenotype'].apply(
                self.extract_omim
            )
        return data

    @staticmethod
    def load_phenotypic_series(filepath: str, colnames: list) -> dict:
        '''Load phenotypic series information from phenotypicSeries.txt
        provided by OMIM.'''
        ps_dict = {}
        with open(filepath, "r") as ps_data:
            for raw_entry in csv.DictReader(
                    (r for r in ps_data if not r.startswith("#")),
                    delimiter="\t",
                    fieldnames=colnames
            ):
                # row with phenotypic series name does not conform to field
                # names
                if raw_entry["ps_number"] not in ps_dict:
                    ps_dict[raw_entry["ps_number"]] = {
                        "name": raw_entry["mim_number"],
                        "syndromes": []
                    }
                else:
                    ps_dict[raw_entry["ps_number"]]["syndromes"].append(
                        {
                            "mim_number": raw_entry["mim_number"],
                            "name": raw_entry["phenotype"]
                        }
                    )
        return ps_dict

    @staticmethod
    def create_phen_omim_to_ps(data):
        '''Create mapping between phenotypic omim ids and phenotypic series.
        '''
        mim_to_ps = collections.defaultdict(list)
        for ps_number, entry in data.items():
            for syndrome in entry["syndromes"]:
                mim_to_ps[syndrome["mim_number"]].append(ps_number)
        return mim_to_ps

    @staticmethod
    def create_phen_to_mim(data):
        '''Create dictionary mapping phenotypic omim ids to list of gene ids.
        '''
        pheno_entries = data.to_dict('records')
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

    @staticmethod
    def load_dataframe(filename, colnames):
        '''Get file from local filesystem.
        If it doesnt exists or cached is false, save from remote url to the
        local path and thereafter read from local filesystem.
        '''
        if not os.path.exists(filename):
            raise RuntimeError(
                "{} does not exist. Add them to the location manually.")
        return pandas.read_table(
            filename, delimiter='\t', comment='#', names=colnames, dtype=str)

    @staticmethod
    def extract_omim(raw_omim):
        '''Extract phenotypic MIM number from the phenotypic string.
        '''
        if pandas.isna(raw_omim):
            return None
        match = RE_OMIM_PHEN.search(raw_omim)
        if match is None:
            return None
        return match.group(1)

    @staticmethod
    def search_table(table, needle, target, default=""):
        '''Search entry in specified column, if not found return default.'''
        try:
            result = table.at[needle, target]
        except KeyError:
            result = default
        return result

    def mim_gene_to_entrez_id(self, mim_gene):
        '''Omim gene id to entrez gene id.'''
        return self.search_table(
            self.indexes["gene_omim"],
            str(mim_gene),
            "entrez_id"
        )

    def entrez_id_to_mim_gene(self, entrez_id):
        '''Entrez gene id to OMIM gene id.'''
        return self.search_table(
            self.indexes["entrez_id"],
            str(entrez_id),
            "mim_number"
        )

    def entrez_id_to_symbol(self, entrez_id):
        '''Entrez gene id to gene symbol (letter code.'''
        return self.search_table(
            self.indexes["entrez_id"],
            str(entrez_id),
            "gene_symbol"
        )

    def mim_pheno_to_mim_gene(self, mim_pheno):
        '''Phenotypic omim id to list of genomic omim ids.'''
        mim_pheno = str(mim_pheno)
        if mim_pheno in self.indexes["pheno_omim"]:
            return [
                mim_gene for mim_gene, _ in
                self.indexes["pheno_omim"][str(mim_pheno)]
            ]
        return []

    def mim_pheno_to_gene(self, mim_pheno):
        '''Convert a phenotype omim id to a gene dictionary object containing
        keys: gene_id, gene_symbol and gene_omim_id.
        '''
        if str(mim_pheno) in self.indexes["pheno_omim"]:
            gene_entry = {}
            for mim_gene, gene_symbol in \
                    self.indexes["pheno_omim"][str(mim_pheno)]:
                entrez_id = self.mim_gene_to_entrez_id(mim_gene)
                gene_entry[mim_gene] = {'gene_id': entrez_id,
                                        'gene_symbol': gene_symbol,
                                        'gene_omim_id': mim_gene}
            return gene_entry

        return {}

    @omim_check
    def omim_id_to_phenotypic_series(self, omim_id: str) -> str:
        '''Translate omim id to phenotypic series id, or if it does
        not exist to the empty string.'''

        ps_label = ""
        # take the first phenotypic series entry or if empty leave the
        # empty string
        index = self.indexes["pheno_series"]
        if omim_id in index and index[omim_id]:
            ps_label = index[omim_id][0]
        return ps_label

    @omim_check
    def _replace_deprecated(self, omim_id: str) -> list:
        '''Replace omim ids that are deprecated or have been moved.'''

        if omim_id in self.files["omim_deprecated_replacement"]:
            return self.files["omim_deprecated_replacement"][omim_id]
        return [omim_id]

    @omim_list
    def replace_deprecated_all(
            self, omim_ids: typing.Union[str, list]
    ) -> list:
        '''Replace all omim ids in the list and deduplicate duplicated
        omim ids.
        Accepts single values and iterables.
        '''
        replaced_ids = [
            str(r)
            for o in omim_ids
            for r in self._replace_deprecated(o)
        ]
        return list(set(replaced_ids))
