'''
Constants
---
Constants used in other scripts. These are mostly interpretations of fields
provided in the Face2Gene jsons.
'''

CHROMOSOMAL_TESTS = [
        'CHROMOSOMAL_MICROARRAY',
        'FISH',
        'KARYOTYPE'
        ]

POSITIVE_RESULTS = [
        'ABNORMAL',
        'ABNORMAL_DIAGNOSTIC',
        'DELETION_DUPLICATION',
        'VARIANTS_DETECTED'
        ]

HGVS_OPS = {
        'SUBSTITUTION': '>',
        'DELETION': 'del',
        'DUPLICATION': 'dup',
        'INSERTION': 'ins',
        'INVERSION': 'inv',
        'DELETION_INSERTION': 'delins',
        'UNKNOWN': ''
        }

HGVS_PREFIX = {
        'CDNA_LEVEL': 'c',
        'PROTEIN_LEVEL': 'p',
        'GENOMIC_DNA_LEVEL': 'g',
        'UNKNOWN': '',
        'RS_NUMBER': ''
        }
