'''
Test case functions
'''
import os
import unittest
from lib.api import phenomizer
from lib.model import json_parser, case
from tests.test_json_loading import BaseMapping


class CaseTest(BaseMapping):
    '''Case method tests.'''

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.phenomizer = phenomizer.PhenomizerService(config=cls.config)

    def load_json(self, name: str) -> dict:
        loaded = json_parser.NewJson.from_file(
            self.get_case_path(name)
        )
        return loaded

    def load_case(self, name: str) -> case.Case:
        case_obj = case.Case(
            self.load_json(name),
            error_fixer=self.error_fixer,
            omim_obj=self.omim,
        )
        case_obj.phenomize(self.phenomizer)
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
                check_status, issues = tcase.check(self.omim)
                if not check_status:
                    print(issues)
                self.assertEqual(check_status, valid, msg)

    def test_multi_diagnosis(self):
        '''Assert that cases with multiple diagnoses are correctly sorted.'''
        tests = [
            ("multi_diagnosis_0.json", True),
            ("multi_diagnosis_1.json", True),
            ("multi_diagnosis_2.json", True),
            ("multi_diagnosis_3.json", True),
            ("multi_diagnosis_4.json", True),
            ("multi_diagnosis_5.json", True),
            ("multi_diagnosis_6.json", True),
            ("multi_diagnosis_7.json", True),
            ("multi_diagnosis_8.json", True),
            ("multi_diagnosis_9.json", True),
            ("multi_diagnosis_10.json", True),
            ("multi_diagnosis_11.json", True),
            ("multi_diagnosis_12.json", True),
            ("multi_diagnosis_13.json", True),
            ("multi_diagnosis_14.json", True),
            ("multi_diagnosis_15.json", True),
            ("multi_diagnosis_16.json", True),
            ("multi_diagnosis_17.json", True),
            ("multi_diagnosis_18.json", True),
            ("multi_diagnosis_19.json", True),
        ]

        for test, result in tests:
            with self.subTest(i=test):
                tcase = self.load_case(test)
                check_stats, issues = tcase.check(self.omim)
                if not check_stats:
                    print(issues)
                self.assertEqual(check_stats, result)
