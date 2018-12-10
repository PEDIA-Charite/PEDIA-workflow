'''
Parse the input args from command line.
'''
import sys
from argparse import ArgumentParser

class PEDIAParser():

    def __init__(
            self
    ):
        super().__init__()
        '''Command line arguments affecting preprocess run behavior.'''
        parser = ArgumentParser(description=(
            "Process f2g provided jsons into a format processable by "
            "classification."))
        parser.add_argument("-s", "--single", help="Process a single json file.")
        parser.add_argument("-l", "--lab", help="Name of the lab you which you define in config.ini.")
        parser.add_argument(
            "--lab-format",
            action='store_true',
            help="Parsing the JSON input with F2G LAB foramt"
        )
        parser.add_argument("-v", "--vcf", help="Path of the real vcf for the case you want to run PEDIA")
        parser.add_argument("--lab-case-id", help="Lab case ID, use it with --lab.")
        parser.add_argument(
            "-o", "--output",
            help="Destination of created old json.",
            default=""
        )
        parser.add_argument(
            "-p", "--pickle",
            help="Start with pickled cases after phenomization."
        )
        parser.add_argument(
            "-e", "--entry",
            help=("Start entrypoint for pickled results. "
                  "Default: pheno - start at phenomization. "
                  "convert - start at old json mapping. "
                  "qc - start at case quality check. "
                  "Used in conjunction with --pickle."),
            default="pheno"
        )
        parser.add_argument(
            "--skip-vcf", action='store_true',
            help="Skip vcf convertion."
        )
        parser.add_argument(
            "-c", "--convert-failed", action='store_true',
            help="Convert cases that failed the first part of quality control."
        )

        parser.add_argument(
            "--aws-format", action='store_true',
            help="Use the AWS format for parsing JSON file"
        )

        self.args = parser.parse_args()

        if self.args.lab and not self.args.lab_case_id:
            parser.print_help()
            sys.exit('Error: No lab case ID! Please provide lab case ID with --case-id')
        if self.args.lab_case_id and not self.args.lab:
            parser.print_help()
            sys.exit('Error: No lab name! Please provide lab name with --lab')
