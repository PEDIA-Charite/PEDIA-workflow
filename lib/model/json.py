import json
import functools
import os
from pprint import pprint

from lib.utils import optional_descent, checkHGVS
from lib.model.hgvs import HGVSModel
from lib.model.syndrome import Syndrome
from lib.constants import CHROMOSOMAL_TESTS, POSITIVE_RESULTS

'''
JSON Models
---
Models to accept json files as provided by Face2Gene. Some conversion and checking facilities will be implemented here.

Class overview:
    JsonFile - Base class handling IO
    OldJson - Old json format for compatibility concerns
    NewJson - Current json format
'''

class JsonFile:
    '''Base json class providing file system operations and basic schema operations.
    This class can be extended to accomodate specific format specifications.
    Extension json classes should implement the following functions for the case class:
    get_case_id - string of case id
    get_variants - list of hgvs variant objects
    get_genome_suggestions - list of syndrome objects
    get_features - list of hpo terms
    get_diagnosis - single syndrome suggestion object
    get_submitter - dict with keys team, email, name
    '''
    def __init__(self, path, corrected_location=None):

        self._correct_dir = corrected_location
        self._path = path
        # our aws bucket is always structured into bucket_dir/cases containing case jsons
        # and bucket_dir/genomic_entries containing mutation information
        base, self._filename = os.path.split(path)
        self._basedir, self._jsondir = os.path.split(base)
        # load a corrected json if it exists and is given
        if self._correct_dir:
            override = os.path.join(self._correct_dir, self._jsondir, self._filename)
            if os.path.exists(override):
                path = override
        self._rawjs = json.load(open(path, 'r'))

    def _load_json(self, directory, entry_id):
        '''Load a json file based on id from specified intermediary directory.
        '''
        filename = '{}.json'.format(entry_id)
        entries_path = os.path.join(self._basedir, directory, filename)
        # override entry path if another corrected path is available
        if self._correct_dir:
            corrected_path = os.path.join(self._correct_dir, directory, filename)
            if os.path.exists(corrected_path):
                entries_path = corrected_path
        if not os.path.exists(entries_path):
            raise OSError("File {} not found".format(genomic_entry))
        return json.load(open(entries_path, 'r'))

    def generate(self):
        '''Get schema of rawjs.'''
        return self._generate_schema(self.raw)

    def _generate_schema(self, data):
        '''Create schema of given raw json data. Some features might not be inferred further.'''
        if isinstance(data, dict):
            return {k:self._generate_schema(v) for k,v in data.items()}
        elif isinstance(data, list):
            res = [self._generate_schema(v) for v in data]
            out = []
            for r in res:
                if isinstance(r, dict):
                    if set(r.keys()) not in [set(o) for o in out]:
                        out.append(r)
                elif isinstance(r, str):
                    if r not in out:
                        out.append(r)
                else:
                    out.append(r)
            return out
        else:
            return ''

    def check(self, false_only=True):
        '''Check schema of current self object and return errors detected.
        '''
        schema =  self.check_schema(self.schema, self.raw)
        self.error = schema
        if false_only:
            return self.filter_schema(schema)
        else:
            return schema

    def check_schema(self, schema, rawjs):
        '''Recursively check the specified schema against the provided schema
        '''
        if isinstance(schema, dict):
            if not isinstance(rawjs,dict):
                return (False, "Not a dict")
            else:
                out = {}
                for k,v in schema.items():
                    if k not in rawjs:
                        out[k] = (False, "No key")
                    else:
                        out[k] = self.check_schema(v,rawjs[k])
                return out
        elif isinstance(schema, list):
            if not isinstance(rawjs,list):
                return (False, "Not a list")
            else:
                if len(schema) == 0:
                    return (True, "")
                out = []
                for v in rawjs:
                    res = []
                    for s in schema:
                        res.append(self.check_schema(s,v))
                    if len(schema) == 1:
                        res = res[0]
                    out.append(res)
                return out
        elif hasattr(schema, '__call__'):
            return schema(rawjs)
        else:
            if rawjs is None:
                return (False, "No value")
            elif schema == '':
                return (True, "")
            elif schema == rawjs:
                return (True, "matches expected string")
            else:
                return (False,"No value")

class OldJson(JsonFile):
    '''Implementing the old Face2Gene json schema.
    '''

    # schema defines how a well formed json version should be defined
    # it can be used to check integrity of the provided json
    schema = {
            'submitter': {
                'team':''
                ,'name':''
            }
            ,'vcf':''
            ,'geneList':[{'gestalt_score':''}]
            ,'features':''
            ,'ranks':''
            ,'genomicData':[
                {   'Mutations':{'HGVS-code':'','Mutation 1':''}
                    ,'Test Information':{'Mutation Type':''}
                }
            ]
    }


    def __init__(self, instdata, omim=None):
        # mapping defines how another json can be translated into the current format
        self.mapping_from_new = {
            'submitter' : {
                'team' : lambda x: x['submitter']['user_team'],
                'name' : lambda x: x['submitter']['user_name']
                },
            # get the key or set the default noVCF
            'vcf' : lambda x: x.get('vcf','noVCF'),
            'geneList' : self._syndrome_to_geneList,
            'features' : lambda x: x['features'],
            'ranks' : lambda x: x['selected_syndromes'],
            'genomicData' : self._genomic_entries_to_genomicData
            }
        self._omim = omim

        if isinstance(instdata, str):
            super().__init__(instdata)
        else:
            super().__init__()
        if isinstance(instdata, NewJson):
            assert omim is not None, 'OMIM object needed for mapping from new json'
            self.map_from_new(instdata)

    def map_from_new(self, newJson):
        mapped = self._map(newJson._rawjs, self.mapping_from_new)
        self._rawjs = mapped
        return mapped

    def _map(self, data, mapping):
        '''Recursively traverse the mapping, but statically reference the originally data.'''
        if isinstance(mapping, dict):
            return {k:self._map(data, v) for k,v in mapping.items()}
        # check that we can call the mapping as a function if it is not another level in our mapping dict
        elif hasattr(mapping, '__call__'):
            return mapping(data)
        else:
            raise TypeError

    def _syndrome_to_geneList(self, data):
        '''Convert a syndrome entry to older geneList format.
        These functions should never need to be used directly.
        '''
        geneList = []
        for s in data['selected_syndromes']:
            # one liner creating omim list, this works via short circuiting in python
            omim_list = not isinstance(s['omim_id'],list) and [s['omim_id']] or s['omim_id']
            omim_genes  = [ self._omim.mim_pheno_to_mim_gene(o) for o in omim_list ]
            omim_genes = set(filter(lambda x: x is not None, omim_genes))
            for go in omim_genes:
                g = {
                        'syndrome_name' : s['syndrome_name'],
                        'combined_score' : s.get('combined_score', None),
                        'feature_score' : s.get('feature_score', None),
                        'gestalt_score' : s.get('gestalt_score', None),
                        'has_mask' : s.get('has_mask', None),
                        'gene_id' : self._omim.mim_gene_to_entrez_id(go),
                        'gene_symbol' : self._omim.mim_gene_to_symbol(go),
                        'gene_omim_id' : go
                        }
                geneList.append(g)
        return geneList

    def _genomic_entries_to_genomicData(self, data):
        '''Convert genomic entries to genomicData format.
        '''

        entry_list = []
        for entry in data['genomic_entries']:
            if entry['test_type'] is None:
                print('Entry is empty.')
                print(entry)
                continue
            hgvs = HGVSModel(entry)
            testInformation = {
                    'Molecular Test' : entry.get('test_type',''),
                    'Notation' : optional_descent(entry, ('variants','variant_information')),
                    'Genotype' : optional_descent(entry, ('variants','zygosity')),
                    'Mutation Type' : optional_descent(entry, ('variants', 'notes')),
                    'Gene Name' : 'gene' in entry and entry['gene'].get('gene_symbol', 'gene_omim_id' in entry['gene'] and 'OMIM-ID:{}'.format(entry['gene']['gene_omim_id']) or '') or ''
                    }

            mutations = {
                    'additional info' : optional_descent(entry, ('variants','notes')),
                    'Build' : optional_descent(entry, ('variants','mutation','chromosome','build')),
                    'result' : entry.get('result', ''),
                    'Inheritance Mode' : optional_descent(entry, ('variants','notes')),
                    'HGVS-code' : hgvs.code,
                    'alternative HGVS-code' : hgvs.alt_code
                    }

            r = {
                    'Test Information' : testInformation,
                    'Mutations' : mutations
                    }
            entry_list.append(r)
        return entry_list

class NewJson(JsonFile):
    '''Implement the new Face2Gene schema as loaded from AWS.
    '''
    schema = {'algo_deploy_version': '',
            'case_id': '',
            'detected_syndromes': [{'combined_score': '',
                'feature_score': '',
                'gestalt_score': '',
                'has_mask': '',
                'omim_id': '',
                'syndrome_name': ''}],
            'documents': [],
            'features': [],
            'genomic_entries': [],
            'selected_syndromes': [{'has_mask': '', 'omim_id': '', 'syndrome_name': ''}],
            'submitter': {'user_email': '', 'user_name': '', 'user_team': ''}}

    def __init__(self, path, *args, **kwargs):
        super().__init__(path, *args, **kwargs)
        # specifiy which fields contain filenames, that have to be loaded afterwards
        self._directive = {
            'genomic_entries':[self._load_genomic]
            }
        # load fields according to directives
        self._load_linked()

    def check(self):
        '''Check whether Json fulfills all provided criteria.
        The criteria are:
            picture has been provided - gestalt_score in detected_syndromes should be greater than 0
            clinical diagnosis - selected_syndromes should not be empty
            single monogenetic disease - not multiple syndromes selected and not multiple pathogenic mutations in different genes
            SNP mutations - no microdel/dup or other large scale aberrations

        Furthermore cases with vcf should not be used for the training of the pipeline.
        '''
        valid = True
        issues = []
        # check maximum gestalt score
        max_gestalt_score = functools.reduce(lambda x,y: max(x,y['gestalt_score']), self._rawjs['detected_syndromes'], 0.0)
        if max_gestalt_score <= 0:
            issues.append('Maximum gestalt score is 0. Probably no image has been provided.')
            valid = False

        # check that only one syndrome has been selected
        if len(self._rawjs['selected_syndromes']) != 1:
            issues.append('{} syndromes have been selected. Only 1 syndrome should be selected for PEDIA inclusion.'.format(len(self._rawjs['selected_syndromes'])))
            valid = False

        # check that molecular information is available at all
        if len(self._rawjs['genomic_entries']) == 0:
            issues.append('No genomic entries available.')
            valid = False

        # check that no structural abnormalities have been detected
        for entry in self._rawjs['genomic_entries']:
            if entry['test_type'] in CHROMOSOMAL_TESTS:
                if entry['result'] in POSITIVE_RESULTS:
                    issues.append('Chromosomal abnormality detected in {} with result {}'.format(entry['test_type'], entry['result']))
                    valid = False

        return valid, issues

    def get_case_id(self):
        return self._rawjs['case_id']

    def get_variants(self):
        models = [ HGVSModel(entry) for entry in self._rawjs['genomic_entries'] ]
        variants = [ v for m in models if m.variants for v in m.variants ]
        return variants

    def get_syndrome_suggestions(self):
        syndromes = [Syndrome(**s) for s in self._rawjs['detected_syndromes'] if s]
        return syndromes

    def get_features(self):
        return self._rawjs['features']

    def get_diagnosis(self):
        diagnosis = [Syndrome(**s,confirmed=True) for s in self._rawjs['selected_syndromes'] if s]
        return diagnosis

    def get_submitter(self):
        submitter = {
                'name' : self._rawjs['submitter']['user_name'],
                'team' : self._rawjs['submitter']['user_team'],
                'email' : self._rawjs['submitter']['user_email']
                }
        return submitter

    def get_vcf(self):
        vcfs = [ d['document_name'] for d in self._rawjs['documents'] if d and d['is_vcf'] ]
        return vcfs

    def _load_linked(self):
        '''
        Call linked with self. properties. This is just a self pointing wrapper around linked.
        '''
        self._rawjs = self._linked(self._rawjs, self._directive)

    def _linked(self, data, load_directive):
        '''Load data according to the provided load directive.
        '''
        if isinstance(data, dict):
            if not isinstance(load_directive, dict):
                raise TypeError
            out = {}
            for k,v in data.items():
                # only load entries, for which we have defined a load directive
                if k in load_directive:
                    out[k] = self._linked(v, load_directive[k])
                else:
                    out[k] = v
            return out
        elif isinstance(data, list):
            if not isinstance(load_directive, list):
                raise TypeError
            # return the original list, if we have no directives defined
            if len(load_directive) == 0:
                return v
            else:
                return [self._linked(v, load_directive[0]) for v in data]
        else:
            if not hasattr(load_directive, '__call__'):
                raise TypeError("Not a loader function.")
            return load_directive(data)

    def _load_genomic(self, genomic_entry):
        return self._load_json('genomics_entries', genomic_entry)
