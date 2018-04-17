'''
Test case functions
'''
import os
import unittest
from lib.model import json, case
from tests.test_json_loading import BaseMapping


class CaseTest(BaseMapping):
    '''Case method tests.'''

    def load_json(self, name: str) -> dict:
        loaded = json.NewJson.from_file(
            self.get_case_path(name)
        )
        return loaded

    def load_case(self, name: str) -> case.Case:
        case_obj = case.Case(
            self.load_json(name),
            error_fixer=self.error_fixer,
            omim_obj=self.omim,
        )
        return case_obj

    def test_check(self):
        tests = [
            (
                "normal.json", True,
                "Normal case did not pass."
            ),
            (
                "diagnosis_many_omim.json", True,
                "Single Syndrome multi omim did not pass."
            ),
            (
                "syndrome_with_card.json", True,
                "Syndrome with available card did not pass."
            )
        ]
        for fname, valid, msg in tests:
            with self.subTest(i=fname):
                tcase = self.load_case(fname)
                check_status, _ = tcase.check(self.omim)
                self.assertEqual(check_status, valid, msg)
