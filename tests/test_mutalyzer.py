'''Mutalyzer API Unittests '''
import os
import unittest
from lib.model import case, json
from lib.api import mutalyzer


INPUT_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data")


class MutalyzerTest(unittest.TestCase):
    '''Testing mutalyzer API calls.'''

    def setUp(self):
        self.mutalyzer = mutalyzer.Mutalyzer()
        input_file = os.path.join(INPUT_PATH, "cases", "51702.json")
        loaded_correct = json.NewJson.from_file(input_file)
        self.cases = [case.Case(loaded_correct)]

    def test_correct_reference_transcripts(self):
        true = ['NM_004992.3:c.473C>T']
        mutalyzer.correct_reference_transcripts(self.cases)
        processed_vars = [str(v) for c in self.cases for v in c.variants]
        self.assertListEqual(
            processed_vars, true,
            "Processed HGVS codes are not equal to reference.")
