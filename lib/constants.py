'''
Constants
---
Constants used in other scripts. These are mostly interpretations of fields
provided in the Face2Gene jsons.
'''

HGVS_ERRORDICT_VERSION = 7

# Bucket name, from where Face2Gene vcf and json files will be downloaded
AWS_BUCKET_NAME = "fdna-pedia-dump"

# tests that count as chromosomal tests, if these are positive, cases will be
# excluded
CHROMOSOMAL_TESTS = [
    'CHROMOSOMAL_MICROARRAY',
    'FISH',
    'KARYOTYPE'
]

# Test result descriptions, that will be counted as positive for our case
# selection criteria
POSITIVE_RESULTS = [
    'ABNORMAL',
    'ABNORMAL_DIAGNOSTIC',
    'DELETION_DUPLICATION',
    'VARIANTS_DETECTED'
]


NEGATIVE_RESULTS = [
    'NORMAL'
    'NORMAL_FEMALE'
    'NORMAL_MALE'
    'NO_SIGNIFICANT_VARIANTS'
]


# Translation of Face2Gene Mutation notation to HGVS operators
HGVS_OPS = {
    'SUBSTITUTION': '>',
    'DELETION': 'del',
    'DUPLICATION': 'dup',
    'INSERTION': 'ins',
    'INVERSION': 'inv',
    'DELETION_INSERTION': 'delins',
    'UNKNOWN': ''
}

# Translation of Description levels in Face2Gene to HGVS sequence types
HGVS_PREFIX = {
    'CDNA_LEVEL': 'c',
    'PROTEIN_LEVEL': 'p',
    'GENOMIC_DNA_LEVEL': 'g',
    'UNKNOWN': '',
    'RS_NUMBER': ''
}

# blacklist HPO illegal hpo terms
ILLEGAL_HPO = [
    'HP:0000006'  # autosomal-dominant inheritance
]
