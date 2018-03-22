'''Mutalyzer API Unittests '''
import os
import unittest

import hgvs.parser

from lib.model import case, json
from lib.api import mutalyzer
from lib.errorfixer import ErrorFixer


INPUT_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data")


class MutalyzerTest(unittest.TestCase):
    '''Testing mutalyzer API calls.'''

    def setUp(self):
        self.mutalyzer = mutalyzer.Mutalyzer()
        input_file = os.path.join(INPUT_PATH, "cases", "normal.json")
        loaded_correct = json.NewJson.from_file(input_file)
        errors = ErrorFixer("", "", version=0)

        self.cases = [case.Case(loaded_correct, errors)]
        self.hgvs = hgvs.parser.Parser()

    def test_correct_reference_transcripts(self):
        true = ['NM_004380.2:c.7302G>A', 'NM_004380.2:c.7302G>A']
        mutalyzer.correct_reference_transcripts(self.cases)
        processed_vars = [str(v) for c in self.cases for v in c.variants]
        self.assertListEqual(
            processed_vars, true,
            "Processed HGVS codes are not equal to reference.")

    def test_rs_to_hgvs(self):
        hgvs_list = self.mutalyzer.get_db_snp_descriptions("rs386834107")
        hgvs_true = [
            'CM000670.2:g.99778780C>T',
            'NC_000008.10:g.100791008C>T',
            'NC_000008.11:g.99778780C>T',
            'NG_007098.2:g.770515C>T',
            'NM_017890.4:c.7603C>T',
            'NM_152564.4:c.7528C>T',
            'NP_060360.3:p.Arg2535Ter',
            'NP_689777.3:p.Arg2510Ter',
            'XP_005250860.1:p.Arg2535Ter',
            'XP_005250858.1:p.Arg2535Ter',
            'XP_005250859.1:p.Arg2534Ter',
            'XP_005250857.1:p.Arg2535Ter',
            'XP_016868598.1:p.Arg2470Ter',
            'XP_016868599.1:p.Arg2130Ter',
            'XP_016868600.1:p.Arg1497Ter',
            'XP_016868601.1:p.Arg1054Ter',
            'XP_011515151.1:p.Arg2509Ter',
            'XP_011515152.1:p.Arg2409Ter',
            'XP_011515153.1:p.Arg1497Ter',
            'XP_011515156.1:p.Arg1128Ter',
            'XP_011515154.1:p.Arg1497Ter',
            'XP_011515155.1:p.Arg2535Ter',
            'XP_011515150.1:p.Arg2534Ter'
        ]
        self.assertListEqual(hgvs_list, hgvs_true)

    def test_check_syntax(self):
        check = self.mutalyzer.check_syntax('NM_001127178.2:c.2005C>T')
        self.assertTrue(check['valid'])

    def test_correct_transcripts(self):
        hgvs_strings = {
            'yolo': 'NM_003002.3:c.274G>T',
            'foker': 'LRG_9t1:c.274G>T',
            'poker': 'chr11:g.111959693G>T',
            'condor': 'NC_000011.9:g.111959693G>T',
            'missing_version': 'NM_001127178:c.2005C>T'
        }
        hgvs_strings_corr = {
            'yolo': 'NM_003002.3:c.274G>T',
            'foker': 'LRG_9t1:c.274G>T',
            'poker': 'chr11:g.111959693G>T',
            'condor': 'NC_000011.9:g.111959693G>T',
            'missing_version': 'NM_001127178.1:c.2005C>T'
        }
        hgvs_strings = {k: [self.hgvs.parse_hgvs_variant(v)]
                        for k, v in hgvs_strings.items()}
        self.mutalyzer.correct_transcripts(hgvs_strings)
        hgvs_strings = {k: str(v[0]) for k, v in hgvs_strings.items()}
        self.assertDictEqual(hgvs_strings, hgvs_strings_corr)
