'''
Test case functions
'''
import os
import unittest
from lib.model import json, case


INPUT_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data")


class CaseTest(unittest.TestCase):
    '''Case method tests.'''

    def setUp(self):
        input_file = os.path.join(INPUT_PATH, "cases", "normal.json")
        loaded_correct = json.NewJson.from_file(input_file)
        self.case = case.Case(loaded_correct)

    def test_check(self):
        self.assertTrue(self.case.check(), "Internal self-check did not pass.")
