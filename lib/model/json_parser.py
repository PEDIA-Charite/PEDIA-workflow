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

from lib.global_singletons import OMIM_INST
from lib.utils import explode_df_column
from lib.vcf_operations import move_vcf
# from lib.utils import optional_descent
from lib.model.hgvs_parser import HGVSModel
from lib import constants


class Directive:
    '''Function to fill a directive'''
    def __init__(self, function: Callable, target: type, intypes: list = []):
        self._target = target
        self._types = intypes
        self._func = self._build_call(function, target)

    def __eq__(self, other):
        if not self._types:
            return True
        return any(isinstance(other, t) for t in self._types)

    @staticmethod
    def _build_call(function, target):
        def call(entry):
            if isinstance(entry, target):
                return entry
            return function(entry)
        return call

    def __call__(self, entry):
        return self._func(entry)


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
    def __init__(
            self,
            data: Union[dict, list],
            path: str = '',
            base_path: str = '',
            override: str = '',
            save_path: str = '',
            file_name: str = '',
            corrected_keys: list = [],
    ):
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
        self._filename = file_name
        self._corrected_keys = corrected_keys

    @classmethod
    def from_file(cls, path: str, corrected_location: str = '') -> 'JsonFile':
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
        override_data = {}
        if corrected_location:
            override = os.path.join(corrected_location, jsondir, filename)
            if os.path.exists(override):
                with open(override, "r") as js_file:
                    override_data = json.load(js_file)

        with open(path, "r") as js_file:
            json_data = json.load(js_file)
            if bool(override_data):
                for key in override_data.keys():
                    json_data[key] = override_data[key]

        LOGGER.debug("Loading json %s", filename)
        # create the parent class
        return cls(
            data=json_data,
            path=path,
            base_path=basedir,
            override=corrected_location,
            corrected_keys=list(override_data.keys()),
        )

    def save_json(
            self,
            save_path: Union[None, str] = None,
            file_name: Union[None, str] = None
    ):
        '''Save the json data contained in self._js
        '''
        save_path = self._save_path if save_path is None else save_path
        file_name = self._filename if file_name is None else file_name

        os.makedirs(save_path, exist_ok=True)
        file_path = os.path.join(save_path, file_name)
        with open(file_path, 'w') as output_json:
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

    def _load_json(self, directory, entry_id, default={}):
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
        if os.path.exists(entries_path):
            with open(entries_path, "r") as entry_file:
                json_data = json.load(entry_file)
        else:
            LOGGER.warning("File %s in %s not found", entry_id, directory)
            json_data = default
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
        if isinstance(load_directive, dict):
            for k, entry in load_directive.items():
                if k not in data:
                    raise IndexError("{} not in data dict.".format(k))
                # do not propagate hidden keys to deeper levels
                data[k] = cls._linked(data[k], entry)
        elif isinstance(load_directive, list):
            data = [
                r for r in
                [
                    cls._linked(v, d)  # hidden keys not used further
                    for d in load_directive
                    for v in data
                ]
                if r
            ]
        elif isinstance(load_directive, Directive):
            if data == load_directive:
                data = load_directive(data)
        return data


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

    def __init__(self, data: dict, save_path: str, file_name: str):
        super().__init__(data, save_path=save_path, file_name=file_name)

    @classmethod
    def from_case_object(cls, case: 'Case', path: str) \
            -> 'OldJson':
        '''Create an old json object from a case entity. This is an alternative
        constructor.
        '''

        LOGGER.debug("Creating OldJson from Case for %s.", case.case_id)
        genomic_data = []
        for model in case.hgvs_models:
            if 'gene_id' in model.gene:
                with open("gene_with_id.txt", "a") as myfile:
                    myfile.write("%s, %s \n" % (case.case_id, model.gene))
            else:
                continue
            data = {
                'Test Information': {
                    'Molecular Test': model.test_type,
                    'Notation': model.variant_info,
                    'Genotype': model.zygosity,
                    'Mutation Type': model.test_type,
                    'Gene Name': model.gene['gene_symbol'],
                    'Gene ID': model.gene['gene_id']
                },
                'Mutations': {
                    'additional info': '',
                    'Build': 'GRCh37',
                    'result': model.result,
                    'Inheritance Mode': '',
                    'HGVS-code': ", ".join(
                        [str(v) for v in model.variants])
                }
            }
            # Remove empty genomic entry
            if data['Test Information']['Gene Name'] != "" \
                    or data['Mutations']['HGVS-code'] != "":
                genomic_data.append(data)

        data = {
            'algo_deploy_version': case.algo_version,
            'case_id': case.case_id,
            'submitter': {
                'user_email': case.submitter['email'],
                'user_team': case.submitter['team'],
                'user_name': case.submitter['name']
            },
            'vcf': case.real_vcf_paths,
            'features': case.features,
            'geneList': case.gene_list,
            'detected_syndromes': case.get_phenomized_list(),
            'genomicData': genomic_data,
            'genomic_entries': case.data.get_js()['genomic_entries'],
            'selected_syndromes': case.data.get_js()['selected_syndromes']
        }
        obj = cls(data, path, "{}.json".format(case.case_id))
        return obj

    def get_case_id(self) -> str:
        '''Return case id of case'''
        return str(self._js["case_id"])


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
        # specifiy which fields contain filenames, that have to be loaded
        # afterwards
        directive = {
            'genomic_entries': [
                Directive(
                    lambda x: self._load_json("genomics_entries", x),
                    target=dict, intypes=[str, int]
                )
            ]
        }
        # load fields according to directives
        self.load_linked(directive)

    def check(self, convert_failed: bool) -> bool:
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
            issues.append(
                ('Maximum gestalt score is 0. Probably no image has '
                 'been provided.')
            )
            valid = False

        # check that no structural abnormalities have been detected
        for entry in self._js['genomic_entries']:
            if entry['test_type'] in constants.CHROMOSOMAL_TESTS:
                if entry['result'] in constants.POSITIVE_RESULTS:
                    issues.append(
                        ('Chromosomal abnormality detected in {} with result '
                         '{}').format(entry['test_type'], entry['result']))
                    valid = False
                    if convert_failed:
                        self._js['genomic_entries'] = []


        return valid, issues

    def get_case_id(self) -> str:
        return str(self._js['case_id'])

    def get_algo_version(self) -> str:
        return str(self._js['algo_deploy_version'])

    def get_js(self):
        return self._js

    def get_genomic_entries(self) -> list:
        return self._js["genomic_entries"]

    def get_variants(self) -> ['HGVSModel']:
        '''Get a list of hgvs objects for variants.
        '''
        models = [HGVSModel(entry)
                  for entry in self._js['genomic_entries']]
        return models

    def get_syndrome_suggestions_and_diagnosis(self) -> pandas.DataFrame:
        '''Return a pandas dataframe containing all suggested syndromes and the
        selected syndroms, which is joined on the table with the confirmed
        column marking the specific entry.
        '''
        # create a dataframe from the list of detected syndromes
        if self._js["detected_syndromes"]:
            syndromes_df = pandas.DataFrame.from_dict(
                self._js['detected_syndromes']
            )
        else:
            syndromes_df = pandas.DataFrame(
                columns=[
                    "omim_id", "gestalt_score", "combined_score",
                    "feature_score", "has_mask", "syndrome_name"
                ]
            )

        # force omim_id to always be a list, required for exploding the df
        syndromes_df['omim_id'] = syndromes_df['omim_id'].apply(
            OMIM_INST.replace_deprecated_all
        )
        # turn omim_list into multiple rows with other properties duplicated
        syndromes_df = explode_df_column(syndromes_df, 'omim_id')
        syndromes_df['omim_id'] = syndromes_df['omim_id'].astype(int)

        # preprare the confirmed diagnosis for joining with the main syndrome
        # dataframe
        if self._js['selected_syndromes']:
            selected_syndromes = [
                dict(
                    s,
                    omim_id=OMIM_INST.replace_deprecated_all(s["omim_id"])
                    or ["0"]
                )
                for s in self._js["selected_syndromes"]
            ]

            selected = pandas.DataFrame.from_dict(selected_syndromes)
            syndromes_df['omim_id'] = syndromes_df['omim_id'].astype(int)

            # create multiple rows from list of omim_id entries duplicating
            # other information
            selected = explode_df_column(selected, 'omim_id')
            selected['omim_id'] = selected['omim_id'].astype(int)
            # add a confirmed diagnosis column
            selected.loc[
                selected["diagnosis"].isin(
                    constants.CONFIRMED_DIAGNOSIS
                ), 'confirmed'
            ] = True
            selected.loc[
                selected["diagnosis"].isin(
                    constants.DIFFERENTIAL_DIAGNOSIS
                ), 'differential'
            ] = True

            # outer join of the syndrome and the confirmed diagnosis
            # pandas.merge has to be used instead of join, because the latter
            # only joins on indices
            syndromes_df = syndromes_df.merge(
                selected, on=['omim_id', 'syndrome_name'], how='outer')
            # set all entries not present in the selected syndromes to not
            # confirmed
            syndromes_df = syndromes_df.fillna({'confirmed': False})
            # merge has_mask
            syndromes_df["has_mask"] = \
                syndromes_df["has_mask_x"].astype(bool) \
                | syndromes_df["has_mask_y"].astype(bool)
            syndromes_df.drop(
                ["has_mask_x", "has_mask_y"], inplace=True, axis=1
            )
            # reset index for continous indexing after the join and explode
            # operations
            syndromes_df = syndromes_df.reset_index(drop=True)
        else:
            # if no syndromes selected, everything is false
            syndromes_df["confirmed"] = False
            syndromes_df["differential"] = False

        syndromes_df['omim_id'] = syndromes_df['omim_id'].astype(int)
        return syndromes_df

    def get_features(self) -> [str]:
        '''Return a list of HPO IDs correponding to entered phenotypic
        features.
        '''
        return [
            h for h in self._js['features'] if h not in constants.ILLEGAL_HPO
        ]

    def get_submitter(self) -> {str: str}:
        '''Return a dictionary containing the submitter name, team and email.
        '''
        submitter = {
            'name': self._js['submitter']['user_name'],
            'team': self._js['submitter']['user_team'],
            'email': self._js['submitter']['user_email']
        }
        return submitter

    def get_vcf(self,
            processed_dir: str = "",
            real_path: str = "") \
            -> [str]:
        '''Get a list of vcf files.
        '''

        if real_path:
            vcfs = [real_path]
        else:
            # vcfs are saved inside documents and marked by is_vcf
            vcfs = [d['document_name']
                    for d in self._js['documents']
                    if d and d['is_vcf']]
            # return empty if no vcfs present
            if not vcfs:
                return []
        case_id = self.get_case_id()
        if os.path.exists(vcfs[0]):
            vcf_path = vcfs[0]
            LOGGER.info("Case %s, VCF file %s is found.",
                        case_id, vcfs[0])
            destination_vcf = os.path.join(processed_dir, case_id + ".vcf.gz")
        else:
            vcf_dir = os.path.join(self._base_dir, "vcfs")
            raw_vcfs = list(os.listdir(vcf_dir))

            # convert and save vcfs to specified location if not already present
            processed_vcfs = [f.strip(".vcf.gz")
                              for f in os.listdir(processed_dir)]
            destination_vcf = os.path.join(processed_dir, case_id + ".vcf.gz")

            case_vcfs = [v for v in raw_vcfs if case_id in v]
            if not case_vcfs:
                LOGGER.info("Case %s, VCF file %s could not be found.",
                            case_id, vcfs[0])
                return []
            vcf = case_vcfs[0]
            vcf_path = os.path.join(vcf_dir, vcf)
        move_vcf(vcf_path, destination_vcf)

        return [destination_vcf]

    def get_detected_syndromes(self) -> [dict]:
        '''Unaltered list of detected syndromes.
        '''
        return self._js['detected_syndromes']

class LabJson(JsonFile):
    '''Implement the new Face2Gene schema as loaded from LAB.
    '''
    # this roughly defines the keys to be expected inside a new format json
    # file
    schema = {
        'lab_case_id': '',
        'case_data': {
            'algo_deploy_version': '',
            'case_id': '',
            'suggested_syndromes': [
                {
                    'syndrome': {
                        'app_valid': '',
                        'omim_id': '',
                        'omim_ids': '',
                        'omim_ps_id': '',
                        'is_group': '',
                        'syndrome_name': ''
                    },
                    'feature_score': '',
                    'gestalt_score': ''
                }
            ],
            'documents': [],
            'features': [],
            'genomic_entries': [],
            'selected_syndromes': [
                {
                    'diagnosis': '',
                    'syndrome': {
                        'app_valid': '',
                        'omim_id': '',
                        'omim_ids': '',
                        'omim_ps_id': '',
                        'is_group': '',
                        'syndrome_name': ''
                    }
                }
            ],
            'posting_user': {'userEmail': '', 'userDisplayName': '', 'userInstitution': ''}
        }
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # specifiy which fields contain filenames, that have to be loaded
        # afterwards
        directive = {
            #'genomic_entries': [
            #    Directive(
            #        lambda x: self._load_json("genomics_entries", x),
            #        target=dict, intypes=[str, int]
            #    )
            #]
        }
        # load fields according to directives
        self.load_linked(directive)

    def check(self, convert_failed: bool) -> bool:
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
        if self._js['case_data']['suggested_syndromes']:
            max_gestalt_score = reduce(
                lambda x, y: max(x, y['gestalt_score']),
                self._js['case_data']['suggested_syndromes'], 0.0)
            if max_gestalt_score <= 0:
                issues.append(
                    ('Maximum gestalt score is 0. Probably no image has '
                     'been provided.')
                )
                valid = False
        else:
            issues.append(
                ('No suggested syndrome. Probably no image has '
                 'been provided.')
            )
            valid = False

        # check that no structural abnormalities have been detected
        if 'genomic_entries' in self._js['case_data']:
            for entry in self._js['case_data']['genomic_entries']:
                if entry['test_type'] in constants.CHROMOSOMAL_TESTS:
                    if entry['result'] in constants.POSITIVE_RESULTS:
                        issues.append(
                            ('Chromosomal abnormality detected in {} with result '
                             '{}').format(entry['test_type'], entry['result']))
                        valid = False
                        if convert_failed:
                            self._js['genomic_entries'] = []


        return valid, issues

    def get_case_id(self) -> str:
        return str(self._js['case_data']['case_id'])

    def get_algo_version(self) -> str:
        return str(self._js['case_data']['algo_version'])

    def get_js(self):
        js = self._js
        js['genomic_entries'] = self.get_genomic_entries()
        js['selected_syndromes'] = self.convert_lab_selected_syndrome(self._js["case_data"]["selected_syndromes"])
        return js

    def get_genomic_entries(self) -> list:
        if 'genomic_entries' in self._js['case_data']:
            return [entry["entry_id"] for entry in self._js['case_data']["genomic_entries"]]
        else:
            return []

    def get_variants(self) -> ['HGVSModel']:
        '''Get a list of hgvs objects for variants.
        '''
        if 'genomic_entries' in self._js['case_data']:
            models = [HGVSModel(entry)
                      for entry in self._js['case_data']['genomic_entries'] if 'variants' in entry]
        else:
            models = []
        return models

    def convert_lab_syndrome(self, syndrome):
        omim_id = syndrome["omim_ids"] if syndrome["is_group"] else syndrome["omim_id"]
        converted = {
                "omim_id": omim_id,
                "syndrome_name": syndrome["syndrome_name"],
                "has_mask": syndrome["app_valid"]
                }
        return converted

    def convert_lab_suggested_syndrome(self, suggested_syndromes):
        converted_syndromes = []
        for syndrome in suggested_syndromes:
            if "syndrome" not in syndrome:
                continue
            converted = self.convert_lab_syndrome(syndrome["syndrome"])
            converted["feature_score"] = syndrome["feature_score"]
            converted["gestalt_score"] = syndrome["gestalt_score"]
            converted["combined_score"] = 0
            converted_syndromes.append(converted)
        return converted_syndromes

    def convert_lab_selected_syndrome(self, selected_syndromes):
        converted_syndromes = []
        for syndrome in selected_syndromes:
            converted = self.convert_lab_syndrome(syndrome["syndrome"])
            converted["diagnosis"] = syndrome["diagnosis"]
            converted_syndromes.append(converted)
        return converted_syndromes

    def get_syndrome_suggestions_and_diagnosis(self) -> pandas.DataFrame:
        '''Return a pandas dataframe containing all suggested syndromes and the
        selected syndroms, which is joined on the table with the confirmed
        column marking the specific entry.
        '''
        # create a dataframe from the list of detected syndromes
        if self._js["case_data"]["suggested_syndromes"]:
            syndromes_df = pandas.DataFrame.from_dict(
                self.convert_lab_suggested_syndrome(self._js["case_data"]['suggested_syndromes'])
            )
        else:
            syndromes_df = pandas.DataFrame(
                columns=[
                    "omim_id", "gestalt_score", "combined_score",
                    "feature_score", "has_mask", "syndrome_name"
                ]
            )

        # force omim_id to always be a list, required for exploding the df
        syndromes_df['omim_id'] = syndromes_df['omim_id'].apply(
            OMIM_INST.replace_deprecated_all
        )
        # turn omim_list into multiple rows with other properties duplicated
        syndromes_df = explode_df_column(syndromes_df, 'omim_id')
        syndromes_df['omim_id'] = syndromes_df['omim_id'].astype(int)

        # preprare the confirmed diagnosis for joining with the main syndrome
        # dataframe
        if self._js['case_data']['selected_syndromes']:
            selected_syndromes = [
                dict(
                    s,
                    omim_id=OMIM_INST.replace_deprecated_all(s["omim_id"])
                    or ["0"]
                )
                for s in self.convert_lab_selected_syndrome(self._js["case_data"]["selected_syndromes"])
            ]

            selected = pandas.DataFrame.from_dict(selected_syndromes)
            syndromes_df['omim_id'] = syndromes_df['omim_id'].astype(int)

            # create multiple rows from list of omim_id entries duplicating
            # other information
            selected = explode_df_column(selected, 'omim_id')
            selected['omim_id'] = selected['omim_id'].astype(int)
            # add a confirmed diagnosis column
            selected.loc[
                selected["diagnosis"].isin(
                    constants.CONFIRMED_DIAGNOSIS
                ), 'confirmed'
            ] = True
            selected.loc[
                selected["diagnosis"].isin(
                    constants.DIFFERENTIAL_DIAGNOSIS
                ), 'differential'
            ] = True

            # outer join of the syndrome and the confirmed diagnosis
            # pandas.merge has to be used instead of join, because the latter
            # only joins on indices
            syndromes_df = syndromes_df.merge(
                selected, on=['omim_id', 'syndrome_name'], how='outer')
            # set all entries not present in the selected syndromes to not
            # confirmed
            syndromes_df = syndromes_df.fillna({'confirmed': False})
            # merge has_mask
            syndromes_df["has_mask"] = \
                syndromes_df["has_mask_x"].astype(bool) \
                | syndromes_df["has_mask_y"].astype(bool)
            syndromes_df.drop(
                ["has_mask_x", "has_mask_y"], inplace=True, axis=1
            )
            # reset index for continous indexing after the join and explode
            # operations
            syndromes_df = syndromes_df.reset_index(drop=True)
        else:
            # if no syndromes selected, everything is false
            syndromes_df["confirmed"] = False
            syndromes_df["differential"] = False

        syndromes_df['omim_id'] = syndromes_df['omim_id'].astype(int)
        return syndromes_df

    def convert_lab_feature(self, features):
        converted = []
        for feature in features:
            if feature['is_present'] == '1':
                converted.append(feature['feature']['hpo_full_id'])
        return converted

    def get_features(self) -> [str]:
        '''Return a list of HPO IDs correponding to entered phenotypic
        features.
        '''
        return [
            h for h in self.convert_lab_feature(self._js['case_data']['selected_features']) if h not in constants.ILLEGAL_HPO
        ]

    def get_submitter(self) -> {str: str}:
        '''Return a dictionary containing the submitter name, team and email.
        '''
        submitter = {
            'name': self._js['case_data']['posting_user']['userDisplayName'],
            'team': self._js['case_data']['posting_user']['userInstitution'],
            'email': self._js['case_data']['posting_user']['userEmail']
        }
        return submitter

    def get_vcf(self,
            processed_dir: str = "",
            real_path: str = "") \
            -> [str]:
        '''Get a list of vcf files.
        '''

        if real_path:
            vcfs = [real_path]
        else:
            if 'documents' not in self._js:
                return []
            # vcfs are saved inside documents and marked by is_vcf
            vcfs = [d['document_name']
                    for d in self._js['documents']
                    if d and d['is_vcf']]
            # return empty if no vcfs present
            if not vcfs:
                return []

        case_id = self.get_case_id()
        if os.path.exists(vcfs[0]):
            vcf_path = vcfs[0]
            LOGGER.info("Case %s, VCF file %s is found.",
                        case_id, vcfs[0])
            destination_vcf = os.path.join(processed_dir, case_id + ".vcf.gz")
        else:
            vcf_dir = os.path.join(self._base_dir, "vcfs")
            raw_vcfs = list(os.listdir(vcf_dir))

            # convert and save vcfs to specified location if not already present
            processed_vcfs = [f.strip(".vcf.gz")
                              for f in os.listdir(processed_dir)]
            destination_vcf = os.path.join(processed_dir, case_id + ".vcf.gz")

            case_vcfs = [v for v in raw_vcfs if case_id in v]
            if not case_vcfs:
                LOGGER.info("Case %s, VCF file %s could not be found.",
                            case_id, vcfs[0])
                return []
            vcf = case_vcfs[0]
            vcf_path = os.path.join(vcf_dir, vcf)
        move_vcf(vcf_path, destination_vcf)

        return [destination_vcf]

    def get_detected_syndromes(self) -> [dict]:
        '''Unaltered list of detected syndromes.
        '''
        return self._js['detected_syndromes']
