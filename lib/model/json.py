'''
JSON Models
---
Models to accept json files as provided by Face2Gene. Some conversion and
checking facilities will be implemented here.

Class overview:
    JsonFile - Base class handling IO
    OldJson - Old json format for compatibility concerns
    NewJson - Current json format
'''

import logging
import json
from typing import Union, Callable
from functools import reduce
import os

import pandas
import filetype

from lib.utils import explode_df_column
from lib.vcf_operations import move_vcf
# from lib.utils import optional_descent
from lib.model.hgvs_parser import HGVSModel
from lib.constants import CHROMOSOMAL_TESTS, POSITIVE_RESULTS


LOGGER = logging.getLogger(__name__)


def reduce_omim(syndrome_dict: dict, f2g: 'Face2Gene') -> dict:
    '''Check if syndrome dict contains a list of omim ids. If yes, query
    Face2Gene Library to see if we can pinpoint a certain omim id.
    '''
    if isinstance(syndrome_dict['omim_id'], list):
        # get a single omim_id from face2gene library
        omim_id = f2g.search_syndrome(
            syndrome_dict['syndrome_name'], syndrome_dict['omim_id'])
        syndrome_dict['omim_id'] = omim_id
    return syndrome_dict


class JsonFile:
    '''Base json class providing file system operations and basic schema
    operations.  This class can be extended to accomodate specific format
    specifications.

    Extension json classes should implement the following functions for the
    case class:
    get_case_id - string of case id
    get_variants - list of hgvs variant objects
    get_genome_suggestions - list of syndrome objects
    get_features - list of hpo terms
    get_diagnosis - single syndrome suggestion object
    get_submitter - dict with keys team, email, name
    '''
    def __init__(self, data: Union[dict, list], path: str='',
                 base_path: str='', override: str='', save_path: str=''):
        '''
        Args:
            data - loaded and deserialized jSON data
            path - Path from which main json has been loaded
            base_path - base folder of original data
            override - base folder for json overrides. this is necessary for
                loading of linked jsons
        '''
        self._js = data
        self._load_path = path
        self._base_dir = base_path
        self._override_dir = override
        self._save_path = save_path

    @classmethod
    def from_file(cls, path: str, corrected_location: str= '') -> 'JsonFile':
        '''Load jSON from file.
        Args:
            path: Path to json file.
            corrected_location: Alternative directory containing file
                                overrides.
        Returns:
            jsonFile object with data loaded in _js.
        '''
        # split path into hierarchy of aws dump
        # our aws bucket is always structured into bucket_dir/cases containing
        # case jsons and bucket_dir/genomic_entries containing mutation
        # information
        # load a corrected json if it exists and is given
        base, filename = os.path.split(path)
        basedir, jsondir = os.path.split(base)
        if corrected_location:
            override = os.path.join(corrected_location, jsondir, filename)
            path = os.path.exists(override) and override or path

        with open(path, "r") as js_file:
            json_data = json.load(js_file)
        # create the parent class
        json_obj = cls(data=json_data, path=path, base_path=basedir,
                       override=corrected_location)

        LOGGER.debug("Loading json %s", filename)
        return json_obj

    def save_json(self):
        '''Save the json data contained in self._js
        '''
        os.makedirs(os.path.split(self._save_path)[0], exist_ok=True)
        with open(self._save_path, 'w') as output_json:
            json.dump(self._js, output_json)

    def generate(self):
        '''Get schema of js.'''
        return self._generate_schema(self.raw)

    @classmethod
    def _generate_schema(cls, data):
        '''Create schema of given raw json data. Some features might not be
        inferred further.'''
        if isinstance(data, dict):
            return {k: cls._generate_schema(v) for k, v in data.items()}
        elif isinstance(data, list):
            partial_result = [cls._generate_schema(v) for v in data]
            out = []
            for part in partial_result:
                if isinstance(part, dict):
                    if set(part.keys()) not in [set(o) for o in out]:
                        out.append(part)
                elif isinstance(part, str):
                    if part not in out:
                        out.append(part)
                else:
                    out.append(part)
            return out
        else:
            return ''

    def check(self, false_only=True):
        '''Check schema of current self object and return errors detected.
        '''
        schema = self.check_schema(self.schema, self.raw)
        self.error = schema
        if false_only:
            return self.filter_schema(schema)
        else:
            return schema

    @classmethod
    def check_schema(cls, schema, data):
        '''Recursively check the specified schema against the provided schema
        '''
        if isinstance(schema, dict):
            if not isinstance(data, dict):
                return (False, "Not a dict")
            else:
                out = {}
                for k, child in schema.items():
                    if k not in data:
                        out[k] = (False, "No key")
                    else:
                        out[k] = cls.check_schema(child, data[k])
                return out
        elif isinstance(schema, list):
            if not isinstance(data, list):
                return (False, "Not a list")
            else:
                if len(schema) == 0:
                    return (True, "")
                out = []
                for entry in data:
                    res = []
                    # multiple entries serve a as an OR option
                    for candidate in schema:
                        res.append(cls.check_schema(candidate, entry))
                    if len(schema) == 1:
                        res = res[0]
                    out.append(res)
                return out
        elif hasattr(schema, '__call__'):
            return schema(data)
        else:
            if data is None:
                return (False, "No value")
            elif schema == '':
                return (True, "")
            elif schema == data:
                return (True, "matches expected string")
            else:
                return (False, "No value")

    def _load_json(self, directory, entry_id):
        '''Load a json file based on id from specified intermediary directory.
        '''
        filename = '{}.json'.format(entry_id)
        entries_path = os.path.join(self._base_dir, directory, filename)
        # override entry path if another corrected path is available
        if self._override_dir:
            corrected_path = os.path.join(
                self._override_dir, directory, filename)
            if os.path.exists(corrected_path):
                entries_path = corrected_path
        if not os.path.exists(entries_path):
            raise OSError("File {} not found".format(entry_id))
        with open(entries_path, "r") as entry_file:
            json_data = json.load(entry_file)
        return json_data

    def load_linked(self, directive:
                    {'str': Callable[[str], Union[dict, list]]}):
        '''Call linked with self. properties. This is just a self pointing
        wrapper around linked.
        Args:
            directive: Define fields, on which loading operations should be
            done
        '''
        self._js = self._linked(self._js, directive)

    @classmethod
    def _linked(cls, data, load_directive):
        '''Load data according to the provided load directive. This will
        recursively traverse the json structure and match it against the
        provided loading directive.
        '''
        if isinstance(data, dict):
            if not isinstance(load_directive, dict):
                raise TypeError
            out = {}
            for k, entry in data.items():
                # only load entries, for which we have defined a load directive
                if k in load_directive:
                    out[k] = cls._linked(entry, load_directive[k])
                else:
                    out[k] = entry
            return out
        elif isinstance(data, list):
            if not isinstance(load_directive, list):
                raise TypeError
            # return the original list, if we have no directives defined
            if len(load_directive) == 0:
                return data
            else:
                return [cls._linked(v, load_directive[0]) for v in data]
        else:
            if not hasattr(load_directive, '__call__'):
                LOGGER.error(data)
                LOGGER.error(load_directive)
                raise TypeError("Not a loader function.")
            return load_directive(data)


class OldJson(JsonFile):
    '''Implementing the old Face2Gene json schema.
    '''

    # schema defines how a well formed json version should be defined
    # it can be used to check integrity of the provided json
    schema = {
        'submitter': {'team': '', 'name': ''},
        'vcf': '',
        'geneList': [{'gestalt_score': ''}],
        'features': '',
        'ranks': '',
        'genomicData': [
            {'Mutations': {'HGVS-code': '', 'Mutation 1': ''},
             'Test Information': {'Mutation Type': ''}
             }
            ]
    }

    def __init__(self, data: dict, save_path: str=''):
        super().__init__(data, save_path=save_path)

    @classmethod
    def from_case_object(cls, case: 'Case', path: str, omim: 'Omim') \
            -> 'OldJson':
        '''Create an old json object from a case entity. This is an alternative
        constructor.
        '''

        LOGGER.debug("Creating OldJson from Case for %s.", case.case_id)
        genomic_data = []
        for model in case.hgvs_models:
            data = {
                'Test Information': {
                    'Molecular Test': model.test_type,
                    'Notation': model.variant_info,
                    'Genotype': model.zygosity,
                    'Mutation Type': model.test_type,
                    'Gene Name': model.gene['gene_symbol']
                },
                'Mutations': {
                    'additional info': '',
                    'Build': 'GRCh37',
                    'result': model.result,
                    'Inheritance Mode': '',
                    'HGVS-code': ", ".join([str(v) for v in model.variants])
                }
            }
            genomic_data.append(data)

        data = {
            'algo_deploy_version': case.algo_version,
            'case_id': case.case_id,
            'submitter': {
                'user_email': case.submitter['email'],
                'user_team': case.submitter['team'],
                'user_name': case.submitter['name']
            },
            'vcf': case.realvcf,
            'features': case.features,
            # maybe disable
            # 'ranks': case.syndromes.to_dict('records'),
            'geneList': case.get_gene_list(omim),
            'detected_syndromes': case.data.get_detected_syndromes(),
            'genomicData': genomic_data
        }
        path = os.path.join(path, '{}.json'.format(case.case_id))
        obj = cls(data, path)
        return obj


class NewJson(JsonFile):
    '''Implement the new Face2Gene schema as loaded from AWS.
    '''
    # this roughly defines the keys to be expected inside a new format json
    # file
    schema = {
        'algo_deploy_version': '',
        'case_id': '',
        'detected_syndromes': [
            {'combined_score': '',
             'feature_score': '',
             'gestalt_score': '',
             'has_mask': '',
             'omim_id': '',
             'syndrome_name': ''}
        ],
        'documents': [],
        'features': [],
        'genomic_entries': [],
        'selected_syndromes': [
            {'has_mask': '', 'omim_id': '', 'syndrome_name': ''}
        ],
        'submitter': {'user_email': '', 'user_name': '', 'user_team': ''}
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # super().from_file(path, *args, **kwargs)
        # specifiy which fields contain filenames, that have to be loaded
        # afterwards
        directive = {
            'genomic_entries': [
                lambda x: self._load_json('genomics_entries', x)
            ]
            }
        # load fields according to directives
        self.load_linked(directive)

    def check(self) -> bool:
        '''Check whether Json fulfills all provided criteria.
        The criteria are:
            picture has been provided - gestalt_score in detected_syndromes
            should be greater than 0
            clinical diagnosis - selected_syndromes should not be empty
            single monogenetic disease - not multiple syndromes selected and
            not multiple pathogenic mutations in different genes
            SNP mutations - no microdel/dup or other large scale aberrations

        Furthermore cases with vcf should not be used for the training of the
        pipeline.
        '''
        valid = True
        issues = []
        # check maximum gestalt score
        max_gestalt_score = reduce(
            lambda x, y: max(x, y['gestalt_score']),
            self._js['detected_syndromes'], 0.0)
        if max_gestalt_score <= 0:
            issues.append('Maximum gestalt score is 0. Probably no image has \
                          been provided.')
            valid = False

        # check that only one syndrome has been selected
        if len(self._js['selected_syndromes']) != 1:
            issues.append(
                '{} syndromes have been selected. Only 1 syndrome should be \
                selected for PEDIA inclusion.'.format(
                    len(self._js['selected_syndromes'])))
            valid = False

        # check that molecular information is available at all
        if len(self._js['genomic_entries']) == 0:
            issues.append('No genomic entries available.')
            valid = False

        # check that no structural abnormalities have been detected
        for entry in self._js['genomic_entries']:
            if entry['test_type'] in CHROMOSOMAL_TESTS:
                if entry['result'] in POSITIVE_RESULTS:
                    issues.append(
                        'Chromosomal abnormality detected in {} with result \
                        {}'.format(entry['test_type'], entry['result']))
                    valid = False

        return valid, issues

    def get_case_id(self) -> str:
        return str(self._js['case_id'])

    def get_algo_version(self) -> str:
        return str(self._js['algo_deploy_version'])

    def get_js(self):
        return self._js

    def get_variants(self, error_fixer: "ErrorFixer") -> ['hgvs']:
        '''Get a list of hgvs objects for variants.
        '''
        models = [HGVSModel(entry, error_fixer)
                  for entry in self._js['genomic_entries']]
        variants = [v for m in models if m.variants for v in m.variants]
        return variants, models

    def get_syndrome_suggestions_and_diagnosis(self) -> pandas.DataFrame:
        '''Return a pandas dataframe containing all suggested syndromes and the
        selected syndroms, which is joined on the table with the confirmed
        column marking the specific entry.
        '''
        # create a dataframe from the list of detected syndromes
        syndromes_df = pandas.DataFrame.from_dict(
            self._js['detected_syndromes'])

        # force omim_id to always be a list, required for exploding the df
        syndromes_df['omim_id'] = syndromes_df['omim_id'].apply(
            lambda x: not isinstance(x, list) and [x] or x)
        # turn omim_list into multiple rows with other properties duplicated
        syndromes_df = explode_df_column(syndromes_df, 'omim_id')
        syndromes_df['omim_id'] = syndromes_df['omim_id'].astype(int)

        # preprare the confirmed diagnosis for joining with the main syndrome
        # dataframe
        selected = pandas.DataFrame.from_dict(
            self._js['selected_syndromes'])

        selected['omim_id'] = selected['omim_id'].apply(
            lambda x: not isinstance(x, list) and [x] or x)
        selected = explode_df_column(selected, 'omim_id')
        # add a confirmed diagnosis column
        selected['confirmed'] = True

        # outer join of the syndrome and the confirmed diagnosis
        # pandas.merge has to be used instead of join, because the latter only
        # joins on indices
        syndromes_df = syndromes_df.merge(
            selected, on=['omim_id', 'syndrome_name'], how='outer')
        # set all entries not present in the selected syndromes to not
        # confirmed
        syndromes_df = syndromes_df.fillna({'confirmed': False})
        # reset index for continous indexing after the join and explode
        # operations
        syndromes_df = syndromes_df.reset_index(drop=True)
        return syndromes_df

    def get_features(self) -> [str]:
        '''Return a list of HPO IDs correponding to entered phenotypic
        features.
        '''
        return self._js['features']

    def get_submitter(self) -> {str: str}:
        '''Return a dictionary containing the submitter name, team and email.
        '''
        submitter = {
            'name': self._js['submitter']['user_name'],
            'team': self._js['submitter']['user_team'],
            'email': self._js['submitter']['user_email']
        }
        return submitter

    def get_vcf(self, processed_dir: str = "data/PEDIA/vcfs/original") \
            -> [str]:
        '''Get a list of vcf files.
        '''
        # vcfs are saved inside documents and marked by is_vcf
        vcfs = [d['document_name']
                for d in self._js['documents']
                if d and d['is_vcf']]
        # return empty if no vcfs present
        if not vcfs:
            return []
        vcf_dir = os.path.join(self._base_dir, "vcfs")
        raw_vcfs = list(os.listdir(vcf_dir))

        # convert and save vcfs to specified location if not already present
        processed_vcfs = [f.strip(".vcf.gz")
                          for f in os.listdir(processed_dir)]
        case_id = self.get_case_id()
        destination_vcf = os.path.join(processed_dir, case_id + ".vcf.gz")

        if case_id not in processed_vcfs:
            case_vcfs = [v for v in raw_vcfs if case_id in v]
            if not case_vcfs:
                LOGGER.warn("Case %s, VCF file %s could not be found.",
                            case_id, vcfs[0])
                return []
            vcf = case_vcfs[0]
            vcf_path = os.path.join(vcf_dir, vcf)
            kind = filetype.guess(vcf_path)
            # get mimetype
            mime = kind.mime if kind is not None else "text"
            move_vcf(vcf_path, destination_vcf, mime)

        return [destination_vcf]

    def get_detected_syndromes(self) -> [dict]:
        '''Unaltered list of detected syndromes.
        '''
        return self._js['detected_syndromes']
