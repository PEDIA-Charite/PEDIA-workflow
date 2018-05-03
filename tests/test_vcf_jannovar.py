import os

from lib import vcf_jannovar
from tests.test_config import BaseConfig


class JannovarTest(BaseConfig):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tmp_path = os.path.join(cls.input_path, "jannovar")
        os.makedirs(cls.tmp_path, exist_ok=True)

    def vcf_path(self, name):
        '''Create vcf output path.'''
        return os.path.join(self.tmp_path, name + ".vcf")

    def test_simple_generation(self):
        '''Generate simple vcf file.'''
        variants = []
        case_id = "123456"
        zygosity = "heterozygous"
        vcf_data = vcf_jannovar.create_vcf(
            variants, zygosity, case_id, self.tmp_path
        )

        vcf_jannovar.write_vcf(vcf_data, self.vcf_path(case_id))
