import os
import shutil
import unittest
from lib import download
from lib.model import config

INPUT_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data")


class DownloadTest(unittest.TestCase):

    def setUp(self):
        self.config = config.ConfigManager(
            os.path.join(INPUT_PATH, "config.ini"))

        self.dl_clean = os.path.join(INPUT_PATH, "test_clean")
        os.makedirs(self.dl_clean)

    def test_download_clean(self):
        download.backup_s3_folder(
            aws_access_key=self.config.aws["access_key"],
            aws_secret_key=self.config.aws["secret_key"],
            download_location=self.dl_clean)

    def tearDown(self):
        pass
        # shutil.rmtree(self.dl_clean)
