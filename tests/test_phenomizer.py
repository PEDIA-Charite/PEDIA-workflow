'''Phenomizer API Unittests'''
import os
import unittest
from lib.model import case, json, config
from lib.api import phenomizer, omim

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

INPUT_PATH = os.path.join(SCRIPTDIR, "data")

CONFIG = config.ConfigManager(os.path.join(INPUT_PATH, "config.ini"))


class PhenomizerTest(unittest.TestCase):
    '''Test phenomization of case files.'''

    def setUp(self):
        input_file = os.path.join(INPUT_PATH, "cases", "51702.json")
        loaded_correct = json.NewJson.from_file(input_file)
        self.case = case.Case(loaded_correct)
        self.phenomizer = phenomizer.PhenomizerService(config=CONFIG)
        self.omim_obj = omim.Omim(
            api_key=CONFIG.omim['api_key'], mimdir=INPUT_PATH,
            use_cached=False)

    def test_phenomize(self):
        self.case.phenomize(self.phenomizer)
        processed = self.case.syndromes
        all_filled = processed.dropna(how="any")
        self.assertGreater(
            max(processed['value_pheno']), 0,
            "Phenomization values are 0 or NaN.")
        self.assertGreater(
            max(processed['value_boqa']), 0,
            "Boqa values are 0 or NaN.")
        self.assertGreater(
            all_filled.shape[0], 0,
            "No intersecting syndromes between F2G data and phenomizer.")

    def test_gene_list_export(self):
        self.case.phenomize(self.phenomizer)
        gene_list_pre = self.case.get_gene_list(
            self.omim_obj, filter_entrez_id=False, recreate=False)
        empty_ones_pre = [g for g in gene_list_pre if g['gene_id'] == '']
        gene_list = self.case.get_gene_list(
            self.omim_obj, filter_entrez_id=True, recreate=True)
        empty_ones = [g for g in gene_list if g['gene_id'] == '']
        self.assertEqual(len(empty_ones), 0,
                         "Entries with gene id remained after filtering.")
        self.assertEqual(
            len(gene_list)+len(empty_ones_pre), len(gene_list_pre),
            "Entry filtering did not correctly remove only empty entries.")

    def tearDown(self):
        # explicitly close session, since requests uses http keepalive
        self.phenomizer.close()
