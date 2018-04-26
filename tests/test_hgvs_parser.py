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
        ]

        for test in tests:
            test_data = self.load_genomic_entry(test)
            model = self.create_model(test_data)
            print(model.variants)


    def test_huge_hgvs(self):
        tests = [
            (
                "huge_hgvs.json",
                ["NM_003220.2:c.826_842delCTGCCTGCAGGGAGACGinsAGGAT"]
            )
        ]

        for test, correct in tests:
            test_data = self.load_genomic_entry(test)
            model = self.create_model(test_data)
            var_strs = [str(v) for v in model.variants]
            print(var_strs)
            self.assertListEqual(var_strs, correct)
