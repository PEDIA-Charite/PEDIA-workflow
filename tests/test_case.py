'''
Test case functions
'''
import os
import unittest
from lib.model import json, case
from tests.test_json_loading import BaseMapping


class CaseTest(BaseMapping):
    '''Case method tests.'''

    def setUp(self):
        pass

    def load_json(self, name: str) -> dict:
        loaded = json.NewJson.from_file(
            self.get_case_path(name)
        )
        return loaded

    def load_case(self, name: str) -> case.Case:
        case_obj = case.Case(
            self.load_json(name),
            error_fixer=self.error_fixer
        )
        return case_obj

    def test_check(self):
        tcase = self.load_case("normal.json")
        check_status, _ = tcase.check(self.omim)
        self.assertTrue(check_status, "Internal self-check did not pass.")

    def test_multi_omim_single_syndrome(self):
        '''Assert that single diagnosis with many omim passes check.'''
        tcase = self.load_case("diagnosis_many_omim.json")
        check_status, check_data = tcase.check(self.omim)
        print(check_data)
        self.assertTrue(
            check_status,
            "Single Diagnosis with multiple OMIM ids did not pass."
        )
