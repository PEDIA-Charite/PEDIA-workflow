'''
Configuration loader to get parameters from a ini file.
---
This is mainly used to define download directories and API access.
'''
import sys
from typing import Union, Iterable, Dict
from configparser import ConfigParser

from lib import errorfixer
from lib.api import mutalyzer, omim, jannovar, phenomizer

from lib.global_singletons import (
    ERRORFIXER_INST, JANNOVAR_INST, OMIM_INST, PHENOMIZER_INST, AWS_INST, LAB_INST
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

        if args.lab:
            LAB_INST.configure(**self.lab_options(args.lab))
        elif not (args.lab or args.single):
            LAB_INST.configure(**self.pedia_lab_options)

        # configure api components
        ERRORFIXER_INST.configure(**self.errorfixer_options)
        JANNOVAR_INST.configure(**self.jannovar_options)
        OMIM_INST.configure(**self.omim_options)
        PHENOMIZER_INST.configure(**self.phenomizer_options)

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
    def pedia_lab_options(self):
        return {
            "lab_id": self["pedia"]["lab_id"],
            "key": self["pedia"]["key"],
            "secret": self["pedia"]["secret"],
        }

    @property
    def download(self):
        return self._download

    def lab_options(self, lab_name):
        if lab_name not in self:
            sys.exit('Error: Lab name is not found in config.ini! Please check if you use the correct lab name in config.ini')
        return {
            "lab_id": self[lab_name]["lab_id"],
            "key": self[lab_name]["key"],
            "secret": self[lab_name]["secret"],
        }

    def parse_input(self, args: "Namespace"):
        download = self["input"].getboolean(
            "download"
        )
        corrected_path = self["input"]["corrected_path"]
        lab_case_id = 0
        lab = ""
        vcf = ""

        if args.lab:
            if args.lab not in self:
                sys.exit('Error: Lab name is not found in config.ini! Please check if you use the correct lab name in config.ini')
            download_path = self[args.lab]["download_path"]
        else:
            download_path = self["pedia"]["download_path"]
        input_files = []

        if args.single:
            download = False
            input_files = [args.single]
        if args.lab and args.case_id:
            lab = args.lab
            lab_case_id = args.case_id

        if args.vcf:
            vcf = args.vcf

        if args.pickle:
            self._picklefiles = args.pickle

        aws_format = True if args.aws_format else False
        return {
            "download": download,
            "corrected_path": corrected_path,
            "download_path": download_path,
            "input_files": input_files,
            "lab_case_id": lab_case_id,
            "vcf": vcf,
            "lab": lab,
            "aws_format": aws_format
        }

    def parse_output(self, args: "Namespace"):
        simulated_vcf_path = self["output"]["simulated_vcf_path"]
        real_vcf_path = self["output"]["real_vcf_path"]
        vcf_config_file = self["output"]["vcf_config_file"]
        converted_path = args.output if args.output else self["output"]["converted_path"]

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
