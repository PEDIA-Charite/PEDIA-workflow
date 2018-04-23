import unittest

from lib.api import face2gene
from tests.test_config import BaseConfig


class Face2GeneTest(BaseConfig):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.f2g = face2gene.F2GLibrary(path=cls.input_path)

    def test_search_syndrome_name(self):
        tests = [
            ("Multiple Congenital Anomalies-Hypotonia-Seizures Syndrome 1; MCAHS1", None),
            ("Multiple Congenital Anomalies-Hypotonia-Seizures Syndrome", None)
        ]
        for test, correct in tests:
            with self.subTest(i=test):
                result = self.f2g.search_syndrome(test)
                print(result)
