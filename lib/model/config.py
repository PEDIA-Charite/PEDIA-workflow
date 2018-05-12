'''
Configuration loader to get parameters from a ini file.
---
This is mainly used to define download directories and API access.
'''
from typing import Union, Iterable, Dict
from configparser import ConfigParser

from lib import errorfixer
from lib.api import mutalyzer, omim, jannovar, phenomizer

from lib.global_singletons import (
    ERRORFIXER_INST, JANNOVAR_INST, OMIM_INST, PHENOMIZER_INST, AWS_INST
)


class PEDIAConfig(ConfigParser):
    '''
    Implement configuration options for objects needed in the PEDIA pipeline.
    '''

    def __init__(
            self,
            args: Union["Namespace", None] = None,
            conf_file: str = "config.ini"
    ):
        self.path = conf_file
        super().__init__()
        self.read(self.path)

        self.logfile_path = self["general"]["logfile"]
        self.dump_intermediate = self["general"].getboolean(
            "dump_intermediate"
        )

        self.input = self.parse_input(args)

        self.output = self.parse_output(args)

        if args.single:
            self.dump_intermediate = False

        # configure api components
        ERRORFIXER_INST.configure(**self.errorfixer_options)
        JANNOVAR_INST.configure(**self.jannovar_options)
        OMIM_INST.configure(**self.omim_options)
        PHENOMIZER_INST.configure(**self.phenomizer_options)
        AWS_INST.configure(**self.aws_options)

    @property
    def errorfixer_options(self):
        return {
            "hgvs_error_file": self["errorfixer"]["error_path"],
            "hgvs_new_errors": self["errorfixer"]["new_error_path"],
            "version": None,
        }

    @property
    def jannovar_options(self):
        return {
            "url": self["jannovar"]["url"],
            "port": int(self["jannovar"]["port"]),
        }

    @property
    def omim_options(self):
        return {
            "mimdir": self["omim"]["mimdir"],
            "mim2gene_hash": self["omim"]["mim2gene_hash"],
            "morbidmap_hash": self["omim"]["morbidmap_hash"],
        }

    @property
    def phenomizer_options(self):
        return {
            "url": self["phenomizer"]["url"],
            "user": self["phenomizer"]["user"],
            "password": self["phenomizer"]["password"],
        }

    @property
    def aws_options(self):
        return {
            "aws_access_key": self["aws"]["access_key"],
            "aws_secret_key": self["aws"]["secret_key"],
        }

    @property
    def download(self):
        return self._download

    def parse_input(self, args: "Namespace"):
        download = self["input"].getboolean(
            "download"
        )
        corrected_path = self["input"]["corrected_path"]
        download_path = self["input"]["download_path"]
        input_files = []

        if args.single:
            download = False
            input_files = [args.single]
        if args.pickle:
            self._picklefiles = args.pickle
        return {
            "download": download,
            "corrected_path": corrected_path,
            "download_path": download_path,
            "input_files": input_files
        }

    def parse_output(self, args: "Namespace"):
        simulated_vcf_path = self["output"]["simulated_vcf_path"]
        real_vcf_path = self["output"]["real_vcf_path"]
        vcf_config_file = self["output"]["vcf_config_file"]
        converted_path = self["output"]["converted_path"]

        valid_case_path = self["output"]["valid_case_path"]
        quality_check_log = self["output"]["quality_check_log"]

        create_log = True

        if args.single:
            create_log = False
        if args.output:
            self._output_base_dir = args.output
        return {
            "simulated_vcf_path": simulated_vcf_path,
            "real_vcf_path": real_vcf_path,
            "vcf_config_file": vcf_config_file,
            "converted_path": converted_path,
            "valid_case_path": valid_case_path,
            "quality_check_log": quality_check_log,
            "create_log": create_log
        }
