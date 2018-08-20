import os

from lib import vcf_jannovar, vcf_operations
from lib.api import jannovar
from tests.test_config import BaseConfig


class JannovarTest(BaseConfig):

    simple_id = "123456"

    tests = [
        (
            [
                "NM_018136.4:c.567_569del",
                "NM_152486.2:c.305+42_305+43insCCCT",
                "XM_005244727.1:c.799C>T",
                "XM_005244727.1:c.964G>C",
            ],
            "heterozygous",
            True,
        ),
        (
            [
                "NM_018136.10:c.567_569del",
                "NM_152486.2:c.305+42_305+43insCCCT",
                "XM_005244727.1:c.799C>T",
                "XM_005244727.1:c.964G>C",
            ],
            "heterozygous",
            False,
        ),
    ]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tmp_path = os.path.join(cls.input_path, "jannovar")
        os.makedirs(cls.tmp_path, exist_ok=True)

        cls.jannovar = jannovar.JannovarClient("localhost", 8888)

    @classmethod
    def vcf_path(cls, name):
        '''Create vcf output path.'''
        return os.path.join(cls.tmp_path, name + ".vcf.gz")

    def test_simple_generation(self):
        '''Generate simple vcf file.'''
        for variants, zygosity, correct in self.tests:
            vcf_data = vcf_jannovar.create_vcf(
                variants, zygosity, self.simple_id, self.tmp_path
            )
            if correct:
                vcf_jannovar.write_vcfdf(
                    vcf_data, self.vcf_path(self.simple_id)
                )

                data = vcf_jannovar.read_vcfdf(self.vcf_path(self.simple_id))

                hgvs_strings = vcf_jannovar.get_hgvs_codes(data)

                for seq in hgvs_strings:
                    self.assertTrue(seq in variants)

                os.remove(self.vcf_path(self.simple_id))
            else:
                self.assertTrue(isinstance(vcf_data, str))

    def test_server_generation(self):
        '''Generate VCF from hgvs strings with Jannovar server.'''
        for variants, zygosity, correct in self.tests:
            vcf_data = self.jannovar.create_vcf(
                variants, zygosity, self.simple_id
            )
            print(vcf_data)

    def test_server_connect(self):
        self.assertTrue(self.jannovar.can_connect())
