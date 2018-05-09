'''Process cases with existing config managing config requirements.'''
import os
import logging
from typing import Union

from lib import errorfixer, qc_logs
from lib.api import omim
from lib.model import case, json_parser, config


LOGGER = logging.getLogger(__name__)


class Processor:
    '''Manage json_parser and case creation.'''

    def __init__(self, config_path: str = "config.ini"):

        self.config = config.ConfigManager(config_path)

        self.error_fixer = errorfixer.ErrorFixer(config=self.config)

        self.omim = omim.Omim(config=self.config)

        self.logs = qc_logs.QCLogs(
            [
                self.config.jsonparser["json_qc_log"],
                self.config.quality["qc_detailed_log"],
            ]
        )

    def get_corrected_dir(self):
        return self.config.preprocess["corrected_location"]

    def get_base_dir(self):
        return self.config.aws["download_location"]

    def get_case_dir(self):
        return os.path.join(self.get_base_dir(), "cases")

    def case_id_to_path(self, case_id: str):
        return os.path.join(self.get_case_dir(), "{}.json".format(case_id))

    def load_json(self, case_id: str) -> Union[json_parser.NewJson, None]:
        case_path = self.case_id_to_path(case_id)
        if os.path.exists(case_path):
            return json_parser.NewJson.from_file(
                case_path, self.get_corrected_dir()
            )
        else:
            LOGGER.info("Could not find file: %s", case_path)
            return None

    def json_to_case(self, obj: json_parser.NewJson) -> case.Case:
        exclude_normal = self.config.preprocess.getboolean(
            "exclude_normal_variants"
        )
        return case.Case(
            obj, self.error_fixer, self.omim, exclude_normal
        )
