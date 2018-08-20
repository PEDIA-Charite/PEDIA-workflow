'''Mutalyzer API Unittests '''
import os
import unittest

import hgvs.parser

from lib.model import case, json_parser
from lib.api import mutalyzer
from lib.errorfixer import ErrorFixer

from tests.test_json_loading import BaseMapping


class MutalyzerTest(BaseMapping):
    '''Testing mutalyzer API calls.'''

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.mutalyzer = mutalyzer.Mutalyzer()
        cls.hgvs = hgvs.parser.Parser()

    def create_case(self, name: str) -> case.Case:
        js_data = json_parser.NewJson.from_file(self.get_case_path(name))
        return case.Case(js_data, self.error_fixer, self.omim)

    def test_correct_reference_transcripts(self):
        test_cases = [
            "normal.json"
        ]
        cases = [self.create_case(n) for n in test_cases]
        true = ['NM_004380.2:c.7302G>A']
        mutalyzer.correct_reference_transcripts(cases)
        processed_vars = [
            str(v) for c in cases for v in c.get_variants()
        ]
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
