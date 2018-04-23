'''
Test loading of jsons from new-style json files.
'''
import os
import unittest
from lib import errorfixer
from lib.model import json_parser, config
from lib.api import mutalyzer, omim, phenomizer, face2gene

from tests.test_config import BaseConfig


ERROR_VERSION = 1


class BaseMapping(BaseConfig):
    '''Base class for mapping related functions.'''

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.error_fixer = errorfixer.ErrorFixer(
            config=cls.config,
            version=ERROR_VERSION
        )

        cls.omim = omim.Omim(config=cls.config)

    def get_case_path(self, name: str) -> str:
        input_file = os.path.join(self.input_path, "cases", name)
        return input_file


class JsonLoadingTest(BaseMapping):
    '''Load different json files.
    Requirements:
        Ensure that genomic entries are correctly loaded.
        Exclusion criteria are correctly applied.
        Data extraction is correct.
    '''

    def setUp(self):
        input_file = os.path.join(self.input_path, "cases", "normal.json")
        error_file = os.path.join(self.input_path, "cases", "error_hgvs.json")
        self.loaded_correct = json_parser.Newjson_parser.from_file(input_file)
        self.loaded_error = json_parser.Newjson_parser.from_file(error_file)

    def test_check(self):
        '''Test whether the file check is passing.'''
        self.assertTrue(
            self.loaded_correct.check(),
            "Checking correct case did not pass."
        )

    def test_get_case_id(self):
        true_id = "1"
        loaded_id = self.loaded_correct.get_case_id()
        self.assertEqual(loaded_id, true_id,
                         "Case ID {} is not equal to expected {}".format(
                             loaded_id, true_id))

    def test_algo_version(self):
        true_algo = "1.0.1"
        loaded_algo = self.loaded_correct.get_algo_version()
        self.assertEqual(loaded_algo, true_algo,
                         "Algo version {} is not equal to expected {}".format(
                             loaded_algo, true_algo))

    def test_get_submitter(self):
        true = {
            'name': 'Anon Ymous',
            'team': 'Üniversity of Testing',
            'email': 'test@email.tt'
        }
        loaded = self.loaded_correct.get_submitter()
        self.assertDictEqual(
            loaded, true,
            "Submitter information {} is not equal to expected {}".format(
                loaded, true))

    def test_get_vcf(self):
        true = []
        loaded = self.loaded_correct.get_vcf()
        self.assertListEqual(true, loaded, "Loaded VCF Info is identical")

    def test_get_variants(self):
        models = self.loaded_correct.get_variants(self.error_fixer)
        variants = [str(v) for m in models for v in m.variants]
        var_correct = ['NM_004380.2:c.7302G>A', 'NM_004380.2:c.7302G>A']
        self.assertListEqual(variants, var_correct)

    def test_get_variants_error(self):
        models = self.loaded_error.get_variants(self.error_fixer)
        variants = [str(v) for m in models for v in m.variants]
        self.assertListEqual(variants, [])
