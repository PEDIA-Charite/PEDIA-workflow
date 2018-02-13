from pprint import pprint
import enum
import os
import csv
import re

import hgvs
import hgvs.parser
import hgvs.validator
import hgvs.exceptions

from lib.utils import list_all_in
from lib.error_solver import ErrorFixer
from lib.api.mutalyzer import Mutalyzer
from lib.constants import HGVS_OPS, HGVS_PREFIX

'''
Parser class to process HGVS Information with diverse content.
'''

RE_PROTEIN = re.compile('\w \((\w+)\)')

def extract_amino(string):
    '''Extract the three letter Code from the Face2Gene Amino acid description.
    They follow the scheme: X (XXX), although nonsense is written as Ter (Stop). So we need to take this into account.
    '''
    m = RE_PROTEIN.search(string)
    m = m and m.group(1) or ''
    # replace stop, which is not a valid three letter code with Ter
    m = m == 'Stop' and 'Ter' or m
    return m

RE_HGVS = re.compile('[gcmnrp]\.\d+', re.IGNORECASE)

MUTALYZER = Mutalyzer()

HGVS_PARSER = hgvs.parser.Parser()
HGVS_VALIDATOR = hgvs.validator.IntrinsicValidator()

ERROR_FIXER = ErrorFixer()

def is_hgvs(string):
    '''Check if a string might be a hgvs code.
    '''
    # remove all whitespace
    string = "".join(string.split())
    match = RE_HGVS.search(string)
    return bool(match)

RE_DOUBLE_BRACKETS = re.compile('([\w_.]+)\([\w_.]+\):([\w_.><]+)\([\w_.>]+\)')

def clean_hgvs(string):
    '''Remove extraneous information from possible hgvs code.
    '''
    # remove any whitespace first
    no_white = "".join(string.split())
    # resolve strings of format NORMAL(PROTEIN):NORMAL(PROTEIN)
    m = RE_DOUBLE_BRACKETS.search(no_white)
    match = m and ":".join([m.group(1), m.group(2)]) or no_white
    return match

def get_multi_field(data, candidate_ids):
    '''Get information from a number of non-empty fields.
    Assert that information in these fields do not differ, as to prevent unforeseen erros further down the line.
    '''
    entry = ''
    for cid in candidate_ids:
        e = cid in data and data[cid] or ''
        if e and entry:
            assert e == entry, '{} and {} are not same'.format(e,entry)
        elif e:
            entry = e
    return entry

class HGVSModel:
    '''Class to model mutation information received from Face2Gene.'''

    def __init__(self, entry_dict, entry_type='new'):
        if entry_type == 'new':
            self._parse_new(entry_dict)
        else:
            raise TypeError('Only parsing of new genomic entry format has been implemented')

    def _parse_new(self, entry_dict):
        '''New gene entry format contains:
        entry_id - entry id of gene entry json file
        test_type
        gene
        result
        variant_type
        variants
        '''
        entry_dict = ERROR_FIXER[entry_dict]
        if not entry_dict:
            print('Skipping empty entry in error dict. It needs to be manually fixed in {} first.'.format(ERROR_FIXER.filename()))
            self.variants = []
            return

        self.entry_id = entry_dict['entry_id']
        self.gene = 'gene' in entry_dict and entry_dict['gene'] or {'gene_id' : '', 'gene_symbol' : '', 'gene_omim_id' : ''}
        self.result = 'result' in entry_dict and entry_dict['result'] or 'UNKNOWN'
        self.test_type = entry_dict['test_type'] or 'UNKNOWN'
        self.variant_type = entry_dict['variant_type'] or 'UNKNOWN'

        variants, candidates = self._parse_variants(entry_dict['variants'])
        # check if valid hgvs codes have been found
        if variants:
            self.variants = variants
        else:
            ERROR_FIXER.add_faulty(entry_dict, candidates)
            self.variants = []

    def _parse_variants(self, variants):
        '''Create variant information from entries in the form of hgvs codes.
        If information is not parseable, create entry in error dictionary which can be filled by hand.

        We initially collect a number of hgvs candidates, which we will try to complete using additional data.

        Params:
            variants Dictionary containing all variant information.

        Returns:
            Tuple containing a list of entered hgvs objects and a list of generated hgvs candidate strings.
        '''
        if not variants:
            return [], []

        hgvs_candidates = []

        variant_information = 'variant_information' in variants and variants['variant_information'] or 'UNKNOWN'

        hgvs_description = 'hgvs_variant_description' in variants and variants['hgvs_variant_description'] or ''
        if is_hgvs(hgvs_description):
            hgvs_candidates.append(hgvs_description)

        notes = 'notes' in variants and variants['notes'] or ''
        if is_hgvs(notes):
            hgvs_candidates.append(notes)

        mutations = self._get_mutations(variants)
        for mutation in mutations:
            hgvs_candidates += self._parse_mutations(mutation, variant_information)
        variants = []
        for cand in hgvs_candidates:
            hg = clean_hgvs(cand)
            if hg:
                try:
                    var = HGVS_PARSER.parse_hgvs_variant(hg)
                    variants.append(var)
                except hgvs.exceptions.HGVSParseError:
                    print('{}: Error parsing {}, uncleaned {}'.format(self.entry_id, hg, cand))
        return variants, hgvs_candidates


    def _get_mutations(self, data):
        '''Get mutation information from mutation fields.
        '''
        if 'mutation' in data:
            mutations = [data['mutation']]
            mutations = [ x for x in mutations if x ]
        elif 'mutation1' in data:
            mutations = [data['mutation1'],data['mutation2']]
        else:
            mutations = {}
        mutations = [ x for x in mutations if x ]
        return mutations

    def _parse_mutations(self, mutation, variant_information):
        '''Get hgvs candidate strings from mutation information.
        '''
        candidates = []
        protein = False
        transcript = 'transcript' in mutation and mutation['transcript'] or ''
        if is_hgvs(transcript):
            candidates.append(transcript)
            # if we a real transcript part, we should be able to split it off
            s = transcript.split(':')[0]
            if not is_hgvs(s):
                transcript = s
            else:
                transcript = ''

        rs_number = 'rs_number' in mutation and mutation['rs_number'] or ''
        if rs_number:
            j = MUTALYZER.getdbSNPDescriptions(rs_number)
            # add the first entry, since we will have a much too large number of entries
            if j:
                candidates.append(j[0])

        # might be something like substitution, deletion etc
        mutation_type = 'mutation_type' in mutation and mutation['mutation_type'] or 'UNKNOWN'

        # get location information
        location_ids = ['location', 'first_amino_position', 'last_amino_position']
        location = 'location' in mutation and mutation['location'] or ''
        # we are given protein codes in format SINGLE_CODE (TRIPLET CODE), since HGVS uses three letter codes
        p1_loc = 'first_amino_position' in mutation and mutation['first_amino_position'] or ''
        p2_loc = 'last_amino_position' in mutation and mutation['last_amino_position'] or ''
        if p1_loc or p2_loc:
            protein = True

        # get entry1 and entry2
        entry1_ids = ['first_amino_acid','original_base']
        entry1 = get_multi_field(mutation, entry1_ids)
        entry2_ids = ['last_amino_acid', 'substituted_base','inserted_bases', 'deleted_bases']
        entry2 = get_multi_field(mutation, entry2_ids)

        if protein:
            hgvs_string = self._build_protein_hgvs(
                    transcript=transcript,
                    prefix=HGVS_PREFIX[variant_information],
                    position1=p1_loc,
                    position2=p2_loc,
                    orig=entry1,
                    operation=HGVS_OPS[mutation_type],
                    sub=entry2)
        else:
            hgvs_string = self._build_hgvs(
                    transcript=transcript,
                    prefix=HGVS_PREFIX[variant_information],
                    position=location,
                    orig=entry1,
                    operation=HGVS_OPS[mutation_type],
                    sub=entry2)

        candidates.append(hgvs_string)

        return candidates

    def _build_hgvs(self, transcript='', prefix='', position='', orig='', operation='', sub=''):
        '''Create hgvs variant string from provided information in these fields.
        This function only fills the provided strings into the fields of interest. It does not try
        to infer any additional information. The created hgvs strings do not have to be correct. They should
        be validated using an additional validator.
        '''
        hgvs_string = "{transcript}:{prefix}.{position}{orig}{operation}{sub}".format(
                transcript=transcript,
                prefix=prefix,
                position=position,
                orig=orig,
                operation=operation,
                sub=sub)
        hgvs_string = "".join(hgvs_string.split())
        return hgvs_string

    def _build_protein_hgvs(self, transcript='', prefix='', position1='', position2='', orig='', operation='', sub=''):
        orig = extract_amino(orig)
        sub = extract_amino(sub)
        if operation == '>':
            operation = ''
            position1 = position1 and position1 or position2
            position2 = position2 and position2 or position1
            assert position1 == position2, 'Position 1 and 2 are not equal'
            position2 = ''
            hgvs_temp = "{transcript}:{prefix}.{orig}{position1}{position2}{sub}"
        else:
            if orig and sub:
                hgvs_temp = "{transcript}:{prefix}.{orig}{position1}_{sub}{position2}{operation}"
            else:
                position1 = position1 and position1 or position2
                position2 = position2 and position2 or position1
                assert position1 == position2, 'Position 1 and 2 are not equal'
                position2 = ''
                hgvs_temp = "{transcript}:{prefix}.{orig}{sub}{position1}{position2}{operation}"
        hgvs_string = hgvs_temp.format(
                transcript=transcript,
                prefix=prefix,
                position1=position1,
                position2=position2,
                orig=orig,
                sub=sub,
                operation=operation)
        hgvs_string = "".join(hgvs_string.split())
        return hgvs_string
