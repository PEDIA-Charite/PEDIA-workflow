import os
import unittest

from lib.api import omim
from lib.model import config

INPUT_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data")


class OmimTest(unittest.TestCase):

    def setUp(self):
        self.config = config.ConfigManager(
            os.path.join(INPUT_PATH, "config.ini"))

    def test_download(self):
        omim_obj = omim.Omim(
            api_key=self.config.omim["api_key"],
            mimdir=INPUT_PATH,
            use_cached=False)
        print(omim_obj.morbidmap)
