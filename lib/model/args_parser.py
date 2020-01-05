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
        '''Command line arguments affecting preprocess run behavior.t '''
        parser = ArgumentParser(description=(
            "Process f2g provided jsons into a format processable by "
            "classification. \n"
            "Example: "
            "python -s tests/data/cases/124.json -v vcfs/124.vcf.gz -o data/output --phenobot-format"))
        parser.add_argument("-s", "--single", help="Path of input json file.")
        parser.add_argument("-v", "--vcf", help="Path of the real vcf for the case you want to run PEDIA")
        parser.add_argument("--vcf-sample-index",
                default=0,
                type=int,
                help="The index of sample in multip vcf. Default: 0")
        parser.add_argument(
            "-o", "--output",
            help="Destination of created old json.",
            default=""
        )

        parser.add_argument(
            "--phenobot-format", action='store_true',
            help="Use the phenobot format for parsing JSON file"
        )
        parser.add_argument(
            "--aws-format", action='store_true',
            help="Use the AWS format for parsing JSON file"
        )
        parser.add_argument("--config", default='config.ini', help="config.ini file.")

        parser.add_argument(
            "--param-c", default=0, type=float,
            help="Parameter C used in PEDIA classifier"
        )
        parser.add_argument(
            "--train_pickle_path", type=str,
            help="Pickle file of training file. The file is in data/train/train_v*.*.p"
        )

        parser.add_argument("-l", "--lab", help="Name of the lab you which you define in config.ini. (LAB api)")
        parser.add_argument("--lab-case-id", help="Lab case ID, use it with --lab. (LAB api)")
        parser.add_argument(
            "--lab-format",
            action='store_true',
            help="Parsing the JSON input with F2G LAB foramt (LAB api)"
        )
        parser.add_argument(
            "--filter-failed", default=False, action='store_true',
            help="Filter cases that failed the first part of quality control."
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
        self.args = parser.parse_args()

        if self.args.lab and not self.args.lab_case_id:
            parser.print_help()
            sys.exit('Error: No lab case ID! Please provide lab case ID with --case-id')
        if self.args.lab_case_id and not self.args.lab:
            parser.print_help()
            sys.exit('Error: No lab name! Please provide lab name with --lab')
