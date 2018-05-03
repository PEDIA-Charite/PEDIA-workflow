import os
import json
import unittest

from lib.model.hgvs_parser import HGVSModel

from tests.test_json_loading import BaseMapping


class HGVSTest(BaseMapping):

    def load_genomic_entry(self, entry_name):
        entry_path = os.path.join(
            self.input_path, "genomics_entries", entry_name
        )

        with open(entry_path, "r") as gfile:
            data = json.load(gfile)
        return data

    def create_model(self, data):
        return HGVSModel(data, error_fixer=self.error_fixer)

    def test_normal_entries(self):
        tests = [
            (
                "deleted_base.json",
                ['NM_000489.4:c.7205delT']
            )

        ]

        for test, correct in tests:
            test_data = self.load_genomic_entry(test)
            model = self.create_model(test_data)
            var_strs = [str(v) for v in model.variants]
            self.assertListEqual(var_strs, correct)

    def test_weird_hgvs(self):
        '''More complicated hgvs variants beyond simple deletions and
        substitutions.'''
        tests = [
            (
                "huge_hgvs.json",
                ["NM_003220.2:c.826_842delCTGCCTGCAGGGAGACGinsAGGAT"]
            ),
            (
                "complex_pos_delins.json",
                ["NM_021828.4(HPSE2):c.1099-4166_1320+840delins23"]
            )
        ]

        for test, correct in tests:
            test_data = self.load_genomic_entry(test)
            model = self.create_model(test_data)
            var_strs = [str(v) for v in model.variants]
            self.assertListEqual(var_strs, correct)

    def test_bracket_entries(self):
        tests = [
            (
                "hgvs_bracket_name.json",
                []
            ),
            (
                "hgvs_bracket_name_reversed.json",
                []
            ),
            (
                "hgvs_double_bracket.json",
                []
            ),
        ]
        for test, correct in tests:
            test_data = self.load_genomic_entry(test)
            model = self.create_model(test_data)
            var_strs = [str(v) for v in model.variants]
            self.assertListEqual(var_strs, correct)
