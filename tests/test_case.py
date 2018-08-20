'''
Test case functions
'''
import os
import unittest
from lib.api import phenomizer, mutalyzer
from lib.model import json_parser, case
from lib.errorfixer import ErrorFixer
from tests.test_json_loading import BaseMapping


class CaseTest(BaseMapping):
    '''Case method tests.'''

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.phenomizer = phenomizer.PhenomizerService(config=cls.config)
        cls.error_fixer = ErrorFixer(
            hgvs_error_file="tests/data/hgvs_errors_real.json",
            hgvs_new_errors="tests/data/hgvs_new_errors.json",
            version=0
        )

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
            ),
            (
                "no_scores_0.json", False,
                "No transmitted detected syndromes"
            ),
            (
                "no_scores_1.json", True,
                "Missing scores because of transfer issues."
            ),
        ]
        for fname, valid, msg in tests:
            with self.subTest(i=fname):
                tcase = self.load_case(fname)
                check_status, issues = tcase.check(self.omim)
                if check_status != valid:
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
            ("multi_diagnosis_11.json", False),  # no syndromes transmitted
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
                if check_stats != result:
                    print(issues)
                self.assertEqual(check_stats, result)
                if issues:
                    self.assertTrue(all(
                        [s["type"] != "MULTI_DIAGNOSIS" for s in issues]
                    ))

    def test_still_multi(self):
        '''Test cases that are still classified as multi diagnosis.'''
        tests = [
            ("still_multi_0.json", True),
            ("still_multi_1.json", True),
            ("still_multi_2.json", True),
            ("still_multi_3.json", True),
            ("still_multi_4.json", True),
            ("still_multi_5.json", True),
            ("still_multi_6.json", True),
            ("still_multi_7.json", True),
            ("still_multi_8.json", True),
            ("still_multi_9.json", True),
            ("still_multi_10.json", True),
            ("still_multi_11.json", True),  # requires corrected hgvs dict
            ("still_multi_12.json", True),
            ("still_multi_13.json", True),
            ("still_multi_14.json", True),
            ("still_multi_15.json", True),
            ("still_multi_16.json", True),
            ("still_multi_17.json", True),
        ]
        for test, result in tests:
            with self.subTest(i=test):
                tcase = self.load_case(test)
                check_stats, issues = tcase.check(self.omim)
                if not check_stats:
                    print(issues)
                self.assertEqual(check_stats, result)

    def test_hpo_no_score(self):
        '''Test cases with hpo features but no scores.'''
        tests = [
            ("no_pheno_0.json", True),
            ("no_pheno_1.json", True),
        ]
        for test, result in tests:
            with self.subTest(i=test):
                tcase = self.load_case(test)
                check_stats, issues = tcase.check(self.omim)
                if not check_stats:
                    print(issues)
                self.assertEqual(check_stats, result)
