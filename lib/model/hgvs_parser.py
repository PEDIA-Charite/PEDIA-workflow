'''
Parser class to process HGVS Information with diverse content.
'''
import re
import logging
from typing import Union

import hgvs
import hgvs.parser
import hgvs.validator
import hgvs.exceptions
import hgvs.assemblymapper
import hgvs.config


from lib.global_singletons import ERRORFIXER_INST, MUTALYZER_INST
from lib.constants import HGVS_OPS, HGVS_PREFIX


LOGGER = logging.getLogger("lib")


# always leave the reference string in hgvs strings
hgvs.config.global_config.formatting.max_ref_length = None


# creation of hgvs objects from hgvs strings
HGVS_PARSER = hgvs.parser.Parser()
# validation of created hgvs objects
HGVS_VALIDATOR = hgvs.validator.IntrinsicValidator()

# REGEX matching proteins denoted by 'X (Triplet code)'
RE_PROTEIN = re.compile(r'\w \((\w+)\)')
# rough matching of hgvs strings
RE_HGVS = re.compile(r'[gcmnrp]\.\d+', re.IGNORECASE)

# match hgvs strings with single brakets containing the gene name
# such as NM_001.1(XYZ1):12del
TRANSCRIPT = r"[A-Za-z]{2,3}_\d+(\.\d+)?"
VARIANT = r"[gcmnrp]\."
POSITION = r"\d+(_\d+)?"

RE_POSITION = re.compile(VARIANT + POSITION)

RE_TRANSCRIPT = re.compile(TRANSCRIPT)

RE_BRACKET = re.compile(r"\(([\w ._]+)\)")

REFSEQ_PROVIDERS = ["NC", "LRG", "NG", "NM", "NR", "NP"]


def is_transcript(cand: str) -> bool:
    '''Check whether string contains a valid transcript number, and given
    reference transcript string is from a valid provider.
    '''
    matched_transcript = bool(RE_TRANSCRIPT.match(cand))
    valid_provider = any(
        [cand.upper().startswith(p) for p in REFSEQ_PROVIDERS]
    )
    return matched_transcript and valid_provider


def is_position(cand: str) -> bool:
    '''Check whether string might describe a HGVS edit with position.'''
    matched_position = bool(RE_POSITION.match(cand))
    return matched_position


def extract_amino(protein_code: str) -> str:
    '''Extract the three letter Code from the Face2Gene Amino acid description.
    They follow the scheme: X (XXX), although nonsense is written as Ter
    (Stop). So we need to take this into account.
    '''
    match = RE_PROTEIN.search(protein_code)
    match = match.group(1) if match else ""
    # replace stop, which is not a valid three letter code with Ter
    return "Ter" if match == "Stop" else match


def is_hgvs(hgvs_candidate: str) -> bool:
    '''Check if a string might be a hgvs code.
    '''
    # remove all whitespace
    hgvs_candidate = "".join(hgvs_candidate.split())
    return bool(RE_HGVS.search(hgvs_candidate))


def clean_hgvs(hgvs_string: str) -> str:
    '''Remove extraneous information from possible hgvs code.
    '''
    # remove any whitespace first
    hgvs_string = "".join(hgvs_string.split())

    # resolve strings with deldel
    hgvs_string = hgvs_string.replace("deldel", "del")
    # exchange wrong "<" for ">"
    hgvs_string = hgvs_string.replace('<', '>')

    # resolve strings of format NORMAL(PROTEIN):NORMAL(PROTEIN)
    # and alternatively less brackets
    if ":" in hgvs_string:
        raw_transcript, raw_posedit = hgvs_string.split(":", 1)
        # iteratively match transcript
        transcripts = []
        while True:
            match = RE_BRACKET.search(raw_transcript)
            if not match:
                if not transcripts or is_transcript(raw_transcript):
                    transcripts.append(raw_transcript)
                break
            if is_transcript(match[1]):
                transcripts.append(match[1])
            raw_transcript = RE_BRACKET.sub("", raw_transcript)

        # iteratively match position
        posedits = []
        while True:
            match = RE_BRACKET.search(raw_posedit)
            if not match:
                if not posedits or is_position(raw_posedit):
                    posedits.append(raw_posedit)
                break
            if is_position(match[1]):
                posedits.append(match[1])

            raw_posedit = RE_BRACKET.sub("", raw_posedit)

        transcript = transcripts[-1]
        posedit = posedits[-1]

        hgvs_string = "{}:{}".format(transcript, posedit)

    return hgvs_string


def get_multi_field(data: dict, candidate_ids: [str]) -> str:
    '''Get information from a number of non-empty fields.
    Assert that information in these fields do not differ, as to prevent
    unforeseen erros further down the line.
    '''
    entry = ''
    for cid in candidate_ids:
        candidate = cid in data and data[cid] or ''
        # check if entry already assigned, check if it matches the previous
        # entry
        if candidate and entry:
            assert candidate == entry, '{} and {} are not same'.format(
                candidate, entry)
        elif candidate:
            entry = candidate
    return entry


def hgvs_identical(
        seqa: hgvs.sequencevariant, seqb: hgvs.sequencevariant
) -> bool:
    '''Compare whether two HGVS strings most probably describe the same
    variant.'''
    # true if all equal
    if str(seqa) == str(seqb):
        return True

    # check if same position has been described with missing info
    acc_eq = seqa.ac == seqb.ac

    posea = seqa.posedit
    poseb = seqb.posedit
    posa = posea.pos
    posb = poseb.pos
    pos_eq = posa == posb

    edita = posea.edit
    editb = poseb.edit
    if edita.type != editb.type:
        edit_eq = False
    elif (hasattr(edita, "ref")
          and (edita.ref and editb.ref)
          and (edita.ref != editb.ref)):
        edit_eq = False
    elif (hasattr(edita, "alt")
          and (edita.alt and editb.alt)
          and (edita.alt != editb.alt)):
        edit_eq = False
    else:
        edit_eq = True
    return acc_eq and pos_eq and edit_eq


class HGVSModel:
    '''Class to model mutation information received from Face2Gene.'''

    def __init__(
            self,
            entry_dict: dict,
    ):
        '''New gene entry format contains:
        entry_id - entry id of gene entry json file
        test_type
        gene
        result
        variant_type
        variants
        See doc/genomic_entry.md for reference on the format of genomic
        entries.
        '''
        self.corrected = False

        self._js = entry_dict
        self.entry_id = entry_dict['entry_id']
        LOGGER.debug('Processing genomic entry %s', self.entry_id)

        self.result = 'result' in entry_dict and entry_dict['result'] or \
            'UNKNOWN'

        gene_top = 'gene' in entry_dict and entry_dict['gene'] or {}
        gene_variant = 'gene' in entry_dict['variants'] \
            and entry_dict['variants']['gene'] or {}
        self.gene = gene_top or gene_variant \
            or {'gene_id': '', 'gene_symbol': '', 'gene_omim_id': ''}
        self.test_type = entry_dict['test_type'] or 'UNKNOWN'
        self.variant_type = entry_dict['variant_type'] or 'UNKNOWN'

        # fix incorrect gene name
        if self.entry_id in ERRORFIXER_INST:
            self._correct_gene_name()

        variants = self._parse_variants(entry_dict['variants'])
        failed = []
        for var in variants:
            checked = MUTALYZER_INST.check_syntax(var)
            if checked and not checked['valid']:
                message = ["{}:{}".format(v['errorcode'], v['message']) for v
                           in checked['messages']['SoapMessage']]
                failed.append(str(var))
                variants.remove(var)
        if failed:
            info = [dict(self._js, message=message)]
            valid_variants = [str(v) for v in variants]
            ERRORFIXER_INST[self.entry_id] = (
                info, valid_variants, failed)
        self.variants = variants

    def get_json(self):
        return self._js

    def _correct_gene_name(self):
        if 'correct_gene' in ERRORFIXER_INST.get_data(self.entry_id):
            self.gene = ERRORFIXER_INST.get_data(
                self.entry_id
            )['correct_gene']

    def _parse_variants(self, variant_dict: dict) -> list:
        '''Create variant information from entries in the form of hgvs codes.
        If information is not parseable, create entry in error dictionary which
        can be filled by hand.

        We initially collect a number of hgvs candidates, which we will try to
        complete using additional data.

        Args:
            variants: Dictionary containing all variant information.

        Returns:
            Tuple containing a list of entered hgvs objects and a list of
            generated hgvs candidate strings.
        '''
        self.zygosity = 'zygosity' in variant_dict \
            and variant_dict['zygosity'] or 'UNKNOWN'

        # might be something like cdna_level etc
        variant_information = 'variant_information' in variant_dict \
            and variant_dict['variant_information'] or 'UNKNOWN'
        self.variant_info = variant_information

        # get information necessary for hgvs assembly
        # this step can be skipped if we already have an override
        if self.entry_id in ERRORFIXER_INST:
            if len(ERRORFIXER_INST[self.entry_id]) > 0:
                variants = ERRORFIXER_INST[self.entry_id]
                variants = [
                    HGVS_PARSER.parse_hgvs_variant(v) for v in variants
                ]
                self.corrected = True
                return variants

        # return empty if variants are empty
        if not variant_dict:
            return []

        # candidates are possible hgvs strings, these are collected from
        # various sources
        hgvs_candidates = []

        hgvs_description = 'hgvs_variant_description' in variant_dict \
            and variant_dict['hgvs_variant_description'] or ''
        if is_hgvs(hgvs_description):
            hgvs_candidates.append(hgvs_description)

        notes = 'notes' in variant_dict and variant_dict['notes'] or ''
        if is_hgvs(notes):
            hgvs_candidates.append(notes)

        mutations = self._get_mutations(variant_dict)
        for mutation in mutations:
            hgvs_candidates += self._parse_mutations(
                mutation, variant_information)

        # try to parse collected hgvs strings
        # only return successfully parsed hgvs strings
        variants = []
        failures = 0
        failed = []
        for candidate in hgvs_candidates:
            cleaned_hgvs = clean_hgvs(candidate)
            if cleaned_hgvs:
                try:
                    var = HGVS_PARSER.parse_hgvs_variant(cleaned_hgvs)
                    if not any([hgvs_identical(v, var) for v in variants]):
                        variants.append(var)
                except hgvs.exceptions.HGVSParseError:
                    failures += 1
                    failed.append(cleaned_hgvs)
        # add failed and partial failues to error dictionaries
        if failures > 0:
            success = [str(v) for v in variants]
            ERRORFIXER_INST[self.entry_id] = ([self._js], success, failed)
        return variants

    def _get_mutations(self, data: dict) -> [dict]:
        '''Get mutation information from mutation fields.
        '''
        if 'mutation' in data:
            mutations = [data['mutation']]
            mutations = [x for x in mutations if x]
        elif 'mutation1' in data:
            mutations = [data['mutation1'], data['mutation2']]
        else:
            mutations = {}
        mutations = [x for x in mutations if x]
        return mutations

    def _parse_mutations(self, mutation: dict, variant_information: str) \
            -> [str]:
        '''Get hgvs candidate strings from mutation information.
        '''
        candidates = []
        protein = False
        transcript = 'transcript' in mutation and mutation['transcript'] or ''
        if is_hgvs(transcript):
            candidates.append(transcript)
            # if we have a real transcript part, we should be able to split it
            # off
            possible = transcript.split(':')[0]
            # if it doesnt contain variant patterns, we use it
            if not is_hgvs(possible):
                transcript = possible
            else:
                transcript = ''

        rs_number = 'rs_number' in mutation and mutation['rs_number'] or ''
        if rs_number:
            j = MUTALYZER_INST.get_db_snp_descriptions(rs_number)
            # add the first entry, since we will have a much too large number
            # of entries
            if j:
                candidates.append(j[0])

        # might be something like substitution, deletion etc
        mutation_type = 'mutation_type' in mutation \
            and mutation['mutation_type'] or 'UNKNOWN'

        # get location information
        location = 'location' in mutation and mutation['location'] or ''
        # we are given protein codes in format SINGLE_CODE (TRIPLET CODE),
        # since HGVS uses three letter codes
        p1_loc = 'first_amino_position' in mutation \
            and mutation['first_amino_position'] or ''
        p2_loc = 'last_amino_position' in mutation \
            and mutation['last_amino_position'] or ''
        if p1_loc or p2_loc:
            protein = True

        # get entry1 and entry2
        entry1_ids = ['first_amino_acid', 'original_base']
        entry1 = get_multi_field(mutation, entry1_ids)
        entry2_ids = [
            'last_amino_acid', 'substituted_base',
            'inserted_bases', 'deleted_bases'
        ]
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

    def _build_hgvs(self,
                    transcript: str='',
                    prefix: str='',
                    position: str='',
                    orig: str='',
                    operation: str='',
                    sub: str='') -> str:
        '''Create hgvs variant string from provided information in these
        fields.  This function only fills the provided strings into the fields
        of interest. It does not try to infer any additional information. The
        created hgvs strings do not have to be correct. They should be
        validated using an additional validator.
        '''
        hgvs_string = \
            "{transcript}:{prefix}.{position}{orig}{operation}{sub}".format(
                transcript=transcript,
                prefix=prefix,
                position=position,
                orig=orig.upper(),
                operation=operation,
                sub=sub.upper())
        hgvs_string = "".join(hgvs_string.split())
        return hgvs_string

    def _build_protein_hgvs(
            self,
            transcript: str='',
            prefix: str='',
            position1: str='',
            position2: str='',
            orig: str='',
            operation: str='',
            sub: str='') -> str:
        '''Create hgvs string for proteins. Since location notation in protein
        is slightly different another tempate has to be used.
        ---
        This function is currently quite bloated. It should be possible to make
        it easier.
        '''
        orig = extract_amino(orig)
        sub = extract_amino(sub)
        if operation == '>':
            operation = ''
            position1 = position1 and position1 or position2
            position2 = position2 and position2 or position1
            assert position1 == position2, 'Position 1 and 2 are not equal'
            position2 = ''
            hgvs_temp = \
                "{transcript}:{prefix}.{orig}{position1}{position2}{sub}"
        else:
            if orig and sub:
                hgvs_temp = ("{transcript}:{prefix}.{orig}"
                             "{position1}_{sub}{position2}{operation}")
            else:
                position1 = position1 and position1 or position2
                position2 = position2 and position2 or position1
                assert position1 == position2, 'Position 1 and 2 are not equal'
                position2 = ''
                hgvs_temp = ("{transcript}:{prefix}.{orig}{sub}{position1}"
                             "{position2}{operation}")
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
