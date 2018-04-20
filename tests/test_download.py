import os
import shutil
import unittest
from lib import download

from tests.test_config import BaseConfig

class DownloadTest(BaseConfig):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        base, data = os.path.split(cls.input_path)
        cls.dl_clean = os.path.join(base, "aws")
        os.makedirs(cls.dl_clean, exist_ok=True)

    def test_download_clean(self):
        download.backup_s3_folder(
            aws_access_key=self.config.aws["access_key"],
            aws_secret_key=self.config.aws["secret_key"],
            download_location=self.dl_clean)

    def tearDown(self):
        pass
        # shutil.rmtree(self.dl_clean)
