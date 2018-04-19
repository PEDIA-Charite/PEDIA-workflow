import os
import unittest
from lib.model import config


class BaseConfig(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.input_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "data"
        )
        cls.config = config.ConfigManager(
            os.path.join(cls.input_path, "config.ini"))
