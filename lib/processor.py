'''Process cases with existing config managing config requirements.'''
import os
import logging
import pickle
from typing import Union

from lib import errorfixer, qc_logs
from lib.visual import progress_bar
from lib.api import omim
from lib.model import case, json_parser, config


LOGGER = logging.getLogger(__name__)


class Processor:
    '''Manage json_parser and case creation.'''

    pickled_attrs = [
        "qc_cases",
        "valid_cases",
        "old_jsons"
    ]

    def __init__(
            self,
            args: Union["Namespace", None] = None,
            config_path: str = "config.ini"
    ):
        self._filepaths = None

        self._jsons = None
        self._qc_jsons = None
        self._valid_jsons = None

        self._cases = None
        self._qc_cases = None
        self._valid_cases = None

        self._old_jsons = None

        self._pickle_intermediate = True

        self.config = config.PEDIAConfig(config_path)

        self.logs = qc_logs.QCLogs(
            [
                self.config.jsonparser["json_qc_log"],
                self.config.quality["qc_detailed_log"],
            ]
        )

        if args:
            self.apply_args(args)

    @property
    def filepaths(self):
        if not self._filepaths:
            self._filepaths = self._obtain_files()
        return self._filepaths

    @property
    def cases(self):
        if not self._cases:
            self._cases = self.create_cases()
        return self._cases

    @property
    def qc_cases(self):
        if not self._qc_cases:
            self._qc_cases = self.check_cases()
        return self._qc_cases

    @property
    def valid_cases(self):
        if not self._valid_cases:
            self._valid_cases = self._create_valid_cases()
        return self._valid_cases

    @property
    def jsons(self):
        if not self._jsons:
            self._jsons = self.create_jsons()
        return self.jsons

    @property
    def qc_jsons(self):
        if not self._qc_jsons:
            self._qc_jsons = self.check_jsons()
        return self._qc_jsons

    @property
    def valid_jsons(self):
        if not self._valid_jsons:
            self._valid_jsons = self._create_valid_jsons()
        return self._valid_jsons()

    @property
    def old_jsons(self):
        if not self._old_jsons:
            self._old_jsons = self._create_old_jsons()
        return self._old_jsons

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
            obj, config=self.config, exclude_benign_variants=exclude_normal
        )

    @progress_bar("Convert and saving old")
    def _create_old_jsons(self):
        for case_obj in self.cases:
            old_json = json_parser.OldJson.from_case_object(
                case_obj,
                self.conversion_path,
                self.config.omim
            )
            # yields into decorator magic
            yield case_obj.case_id, old_json

    def create_old_json(self):
        return self._create_old_jsons()

    def save_old_jsons(self, old_jsons, path):
        for old_json in old_jsons:
            old_json.save_json(save_path=path)

    def save_all_old_jsons(self):
        self.save_old_jsons(
            [j for _, j in self.qc_jsons], self.conversion_path
        )

    def save_valid_old_jsons(self):
        passing_keys = [c.case_id for c in self.valid_cases]
        valid_old = [
            old_json
            for case_id, old_json in self.old_jsons
            if case_id in passing_keys
        ]
        self.save_old_jsons(valid_old, self.valid_conversion_path)

    @progress_bar("Creating vcf files")
    def _create_vcf(self):
        for case_obj in self.cases:
            simulated = self.simulated_vcf_paths
            real = self.real_vcf_paths
            # yields into decorator magic
            yield {
                "case_id": case_obj.case_id,
                "simulated": simulated,
                "real": real
            }

    def create_vcf_config(self):
        all_vcf = self._create_vcf()
        return all_vcf

    @progress_bar("QC cases")
    def _check_cases(self):
        for case_obj in self.cases:
            valid = case_obj.check(self.config.omim)
            yield valid, case_obj

    def _create_valid_cases(self):
        return [c for v, c in self.qc_cases if v[0]]

    @progress_bar("QC jsons")
    def _check_jsons(self):
        for json_obj in self.jsons:
            valid = json_obj.check()
            yield valid, json_obj

    def _create_valid_jsons(self):
        return [j for v, j in self.qc_jsons if v[0]]

    def _load_pickles(self, pickle_arg: str):
        paths = [p.strip() for p in pickle_arg.split(",")]
        for pickle_path in paths:
            # get filename without extension
            attrname = os.path.splitext(
                os.path.split(pickle_path)[0]
            )[0]
            setattr(self, "_"+attrname, pickle.load(pickle_path))

    def apply_args(self, args: "Namespace"):
        '''Apply arguments to files processing.'''
        if args.single:
            # turn on single file processing
            self._filepaths = [args.single]
            self._pickle_intermediate = False

        if args.output:
            self._conversion_path = [args.output]

        if args.pickle:
            # load pickled data into specific fields
            self._load_pickles(args.pickle)

    def save_pickles(self, path: str = "."):
        '''Save pickled attributes to specified path.'''
        if not self._pickle_intermediate:
            LOGGER.debug("Pickling disabled. Skipping.")
            return

        for attr in self.pickled_attrs:
            obj = getattr(self, "_"+attr)
            if not obj:
                LOGGER.warning("Object %s empty. Not pickling.", attr)
                continue
            pickle_path = os.path.join(path, attr+".p")
            with open(pickle_path, "wb") as pfile:
                pickle.dump(obj, pfile)
