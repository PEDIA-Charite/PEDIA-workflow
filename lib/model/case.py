'''
Case model created from json files.
'''
import logging
from typing import Union

import csv
import subprocess
import tempfile
import os

import pandas

from lib.model.json_parser import OldJson, NewJson
from lib.vcf_operations import move_vcf
from lib import constants


LOGGER = logging.getLogger(__name__)


class Case:
    '''
    Exposes the following properties:
    case_id - Unique identifier for the case in Face2Gene
    variants - list of hgvs objects describing valid hgvs variants
    analyzed_syndromes - dictionary of syndromes using omim id as key
    features - list of hpo terms
    diagnosis - list of syndromes selected as diagnosis
    submitter - submitter information containing fields for email, name, team
    realvcf - list of vcf filenames
    '''

    def __init__(
            self, data: Union[OldJson, NewJson],
            error_fixer: "ErrorFixer",
            omim_obj: "Omim",
            exclude_benign_variants: bool = True
    ):
        self.algo_version = data.get_algo_version()
        self.case_id = data.get_case_id()

        # get both the list of hgvs variants and the hgvs models used in the
        # parsing
        self.hgvs_models = data.get_variants(error_fixer)
        self.syndromes = data.get_syndrome_suggestions_and_diagnosis(omim_obj)
        self.features = data.get_features()
        self.submitter = data.get_submitter()
        self.realvcf = data.get_vcf()
        self.gene_scores = None
        # also save the json object to easier extract information from the
        # new format
        self.data = data
        LOGGER.debug("Creating case %s", self.case_id)

        # query settings
        self.exclude_benign_variants = exclude_benign_variants

    def phenomize(self, pheno: 'PhenomizerService') -> bool:
        '''Add phenomization information to genes from boqa and phenomizer.
        Args:
            omim: Omim object to handle id translation.
            pheno: PhenomizerService to handle API calls for phenomizer and
                   boqa.
        '''
        pheno_boqa = pheno.disease_boqa_phenomize(self.get_features())

        pheno_boqa.index = pheno_boqa.index.astype(int)
        # merge pheno and boqa scores dataframe with our current syndromes
        # dataframe which contains face2gene scores
        self.syndromes["omim_id"] = self.syndromes["omim_id"].astype(int)
        self.syndromes = self.syndromes.merge(
            pheno_boqa, left_on='omim_id', how='outer', right_index=True)
        self.syndromes.reset_index(drop=True, inplace=True)

        self.syndromes.rename(
            columns={
                "value_pheno": "pheno_score",
                "value_boqa": "boqa_score",
            },
            inplace=True
        )

        # fill nans created by merge
        self.syndromes.fillna(
            {
                'combined_score': 0.0,
                'feature_score': 0.0,
                'gestalt_score': 0.0,
                'pheno_score': 0.0,
                'boqa_score': 0.0,
                'syndrome_name': '',
                'confirmed': False,
                'has_mask': False,
                'gene-symbol': '',
                'gene-id': '',
                'disease-name_boqa': '',
                'disease-name_pheno': '',
                'disease-id_boqa': '',
                'disease-id_pheno': '',
            },
            inplace=True
        )

        LOGGER.debug("Phenomization case %s success", self.case_id)

        return True

    def get_gene_list(self, omim: "Omim") -> [dict]:
        '''Get a list of genes from the detected syndrome by inferring
        gene phenotype mappings from the phenomizer and OMIM.
        '''
        syndromes = self.syndromes.to_dict("records")

        phenotypic_series_mapping = {}
        for syndrome in syndromes:

            disease_id = syndrome["omim_id"]

            phenotypic_series = omim.omim_id_to_phenotypic_series(
                str(disease_id)
            ) or str(disease_id)

            syndrome_name = (
                syndrome["syndrome_name"]
                or syndrome["disease-name_pheno"]
                or syndrome["disease-name_boqa"]
            )

            genes = list(omim.mim_pheno_to_gene(disease_id).values())
            if syndrome["gene-id"]:
                genes += [
                    {
                        "gene_id": eid,
                        "gene_symbol": omim.entrez_id_to_symbol(eid),
                        "gene_omim_id": omim.entrez_id_to_mim_gene(eid)
                    }
                    for eid in syndrome["gene-id"].split(", ")
                    if eid not in [g["gene_id"] for g in genes]
                ]

            for gene in genes:
                if not gene["gene_id"]:
                    continue

                # uniqueness constraint on phenotypic series and
                # gene_id
                key = "{}|{}".format(phenotypic_series, gene["gene_id"])
                update_data = dict(
                    {
                        "disease_id": disease_id,
                        "phenotypic_series": phenotypic_series,
                        "syndrome_name": syndrome_name,
                        "gestalt_score": syndrome["gestalt_score"],
                        "feature_score": syndrome["feature_score"],
                        "combined_score": syndrome["combined_score"],
                        "pheno_score": syndrome["pheno_score"],
                        "boqa_score": syndrome["boqa_score"]
                    },
                    **gene
                )
                if key in phenotypic_series_mapping:
                    # use the largest scores of two identical mappings
                    phenotypic_series_mapping[key] = {
                        k: (max(v, update_data[k])
                            if not isinstance(v, str)
                            else update_data[k] or v)
                        for k, v in phenotypic_series_mapping[key].items()
                    }
                else:
                    phenotypic_series_mapping[key] = update_data

        return list(phenotypic_series_mapping.values())

    def pathogenic_gene_in_gene_list(
            self,
            omim: Union[None, "Omim"] = None
    ) -> (bool, list):
        '''Check whether diagnosed genetic mutation is pathogenic gene.'''
        variant_gene_names = [
            v.gene["gene_id"] for v in self.get_hgvs_models()
        ]
        gene_list = self.get_gene_list(omim=omim)
        gene_list_ids = [g["gene_id"] for g in gene_list]
        status = [
            (v in gene_list_ids, v)
            for v in variant_gene_names
        ]
        return status

    def check(self, omim: Union[None, "Omim"] = None) -> bool:
        '''Check whether Case fulfills all provided criteria.

        The criteria are:
            picture has been provided - gestalt_score in detected_syndromes
            should be greater than 0
            clinical diagnosis - selected_syndromes should not be empty
            single monogenetic disease - not multiple syndromes selected and
            not multiple pathogenic mutations in different genes
            SNP mutations - no microdel/dup or other large scale aberrations

        Note that this function is more definitive than the json level check,
        as the validity of hgvs parsing has already been established.
        '''
        valid = True
        issues = []
        scores = [
            "gestalt_score", "feature_score", "pheno_score", "boqa_score"
        ]
        max_scores = {s: max(self.syndromes[s]) for s in scores}
        zero_scores = [s for s, n in max_scores.items() if n <= 0]
        # check maximum gestalt score
        if zero_scores:
            issues.append(
                {
                    "type": "MISSING_SCORES",
                    "data": max_scores
                }
            )
            valid = False

        # check that only one syndrome has been selected
        diagnosis = self.syndromes.loc[self.syndromes['confirmed']]
        if len(diagnosis) < 1:
            issues.append(
                {
                    "type": "NO_DIAGNOSIS",
                }
            )
            valid = False

        # check whether multiple diagnoses are in same phenotypic series
        if omim:
            diagnosis_series = [
                omim.omim_id_to_phenotypic_series(str(d)) or str(d)
                for d in diagnosis["omim_id"]
            ]
            if len(set(diagnosis_series)) > 1:
                issues.append(
                    {
                        "type": "MULTI_DIAGNOSIS",
                        "data": {
                            "orig": list(diagnosis["omim_id"]),
                            "series": diagnosis_series
                        }
                    }
                )
                valid = False
        else:
            LOGGER.warning("No omim object. Some checks will not run.")

        if not self.get_variants():
            raw_entries = self.data.get_genomic_entries()
            if not len(raw_entries):
                issues.append(
                    {
                        "type": "NO_GENOMIC",
                    }
                )
            else:
                issues.append(
                    {
                        "type": "MALFORMED_HGVS",
                        "data": raw_entries
                    }
                )
            valid = False

        return valid, issues

    def check_vcf(self):
        issues = []
        valid = True
        # check simulated vcf is correct
        if hasattr(self, "vcf"):
            if isinstance(self.vcf, str):
                issues.append(
                    {
                        "type": "VCF_ERROR",
                        "data": self.vcf
                    }
                )
                valid = False
        else:
            valid = False
        return valid, issues

    def get_variants(
            self,
            exclusion: Union[bool, None] = None
    ) -> ["hgvs"]:
        '''Get list of variants from all hgvs models.
        Params:
            exclusion - Genomic entries marked explicitly as
            normal are excluded from the returned list.
        '''
        variants = [
            v for m in self.get_hgvs_models(exclusion=exclusion)
            for v in m.variants
        ]
        return variants

    def get_benign_excluded(self) -> int:
        '''Get number of variants excluded by benign filters.'''
        return len(self.get_variants(exclusion=False)) \
            - len(self.get_variants(exclusion=True))

    def get_hgvs_models(
            self,
            exclusion: Union[bool, None] = None
    ) -> ["HGVSModel"]:
        '''Get list of hgvs models, containing all processed information
        on hgvs variants.'''
        exclusion = exclusion if exclusion is not None \
            else self.exclude_benign_variants
        # return all models if no exclusion parameter
        if not exclusion:
            return self.hgvs_models

        models = []
        for model in self.hgvs_models:
            if model.result in constants.NEGATIVE_RESULTS:
                LOGGER.debug(
                    ("NEGATIVE_RESULT Case %s Genomic entry %s marked as "
                     "normal and will be ignored"),
                    self.case_id, model.entry_id)
            else:
                models.append(model)
        return models

    def get_features(self):
        '''
        Get list of features. Exclude illegal HPO terms for the phenomizer.
        '''
        return [
            h for h in self.features
            if h not in constants.ILLEGAL_HPO
        ]

    def get_syndrome_list(self):
        '''Get list of syndromes from syndrome table'''
        syndrome_list = self.syndromes.to_dict("records")
        return syndrome_list

    def eligible_training(self) -> bool:
        '''Eligibility of case for training.
        Exclusion criteria are:
            available real vcf
        '''
        return not self.realvcf

    def get_vcf(self) -> list:
        '''Return vcf files. Does not include generated vcf files.'''
        return self.realvcf

    def create_vcf(self, path: str) -> pandas.DataFrame:
        '''Generates vcf dataframe. If an error occurs the error message is returned.
        '''
        with tempfile.NamedTemporaryFile(mode="w+", dir=path) as hgvsfile:
            for v in self.get_variants():
                hgvsfile.write(str(v) + "\n")
            hgvsfile.seek(0)
            with tempfile.NamedTemporaryFile(mode="w+", dir=path, suffix=".vcf") as vcffile:
                try:
                    process = subprocess.run(["java", "-jar", 'data/jannovar/jannovar_0.25/jannovar-cli-0.25-SNAPSHOT.jar', "hgvs-to-vcf", "-d",
                                              'data/jannovar/jannovar_0.25/data/hg19_refseq.ser', "-i", hgvsfile.name, "-o", vcffile.name, "-r", "data/referenceGenome/data/human_g1k_v37.fasta"], check=True, universal_newlines=True, stderr=subprocess.PIPE)
                except subprocess.CalledProcessError as e:
                    return str(e)
                columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                           'QUAL', 'FILTER', 'INFO', 'FORMAT', self.case_id]
                df = pandas.read_table(
                    vcffile.name, sep='\t', comment='#', names=columns)
                df.ALT.fillna("NA", inplace= True)
                if any(df.ALT == '<ERROR>'):
                    return(process.stderr)
                if self.hgvs_models[0].zygosity.lower() == 'hemizygous':
                    genotype = '1'
                elif self.hgvs_models[0].zygosity.lower() == 'homozygous':
                    genotype = '1/1'
                elif self.hgvs_models[0].zygosity.lower() == 'heterozygous' or self.hgvs_models[0].zygosity.lower() == 'compound heterozygous':
                    genotype = '0/1'
                else:
                    genotype = '0/1'
                df[self.case_id] = genotype
                df['FORMAT'] = 'GT'
                df['INFO'] = ['HGVS="' + str(v) + '"' for v in self.get_variants()]
                df = df.sort_values(by=['#CHROM', "POS"])
                df = df.drop_duplicates()
                return df

    def dump_vcf(self, path: str, recreate: bool = False) -> None:
        '''Dumps vcf file to given path. Initializes vcf generation if none has yet been created.
        Created vcf is saved to self.vcf.
        '''
        if hasattr(self, 'vcf') and not recreate:
            if isinstance(self.vcf, str):
                LOGGER.debug(
                    "VCF generation for case %s failed. Error message:%s", self.case_id, self.vcf)
            else:
                outputpath = os.path.join(path, self.case_id + '.vcf')
                # add header to vcf
                with open(outputpath, 'w') as outfile:
                    outfile.write(
                        '##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    outfile.write(
                        '##contig=<ID=1,assembly=b37,length=249250621>\n')
                    outfile.write(
                        '##contig=<ID=2,assembly=b37,length=243199373>\n')
                    outfile.write(
                        '##contig=<ID=3,assembly=b37,length=198022430>\n')
                    outfile.write(
                        '##contig=<ID=4,assembly=b37,length=191154276>\n')
                    outfile.write(
                        '##contig=<ID=5,assembly=b37,length=180915260>\n')
                    outfile.write(
                        '##contig=<ID=6,assembly=b37,length=171115067>\n')
                    outfile.write(
                        '##contig=<ID=7,assembly=b37,length=159138663>\n')
                    outfile.write(
                        '##contig=<ID=8,assembly=b37,length=146364022>\n')
                    outfile.write(
                        '##contig=<ID=9,assembly=b37,length=141213431>\n')
                    outfile.write(
                        '##contig=<ID=10,assembly=b37,length=135534747>\n')
                    outfile.write(
                        '##contig=<ID=11,assembly=b37,length=135006516>\n')
                    outfile.write(
                        '##contig=<ID=12,assembly=b37,length=133851895>\n')
                    outfile.write(
                        '##contig=<ID=13,assembly=b37,length=115169878>\n')
                    outfile.write(
                        '##contig=<ID=14,assembly=b37,length=107349540>\n')
                    outfile.write(
                        '##contig=<ID=15,assembly=b37,length=102531392>\n')
                    outfile.write(
                        '##contig=<ID=16,assembly=b37,length=90354753>\n')
                    outfile.write(
                        '##contig=<ID=17,assembly=b37,length=81195210>\n')
                    outfile.write(
                        '##contig=<ID=18,assembly=b37,length=78077248>\n')
                    outfile.write(
                        '##contig=<ID=19,assembly=b37,length=59128983>\n')
                    outfile.write(
                        '##contig=<ID=20,assembly=b37,length=63025520>\n')
                    outfile.write(
                        '##contig=<ID=21,assembly=b37,length=48129895>\n')
                    outfile.write(
                        '##contig=<ID=22,assembly=b37,length=51304566>\n')
                    outfile.write(
                        '##contig=<ID=X,assembly=b37,length=155270560>\n')
                    outfile.write(
                        '##contig=<ID=Y,assembly=b37,length=59373566>\n')

                self.vcf.to_csv(outputpath, mode='a', sep='\t', index=False,
                                header=True, quoting=csv.QUOTE_NONE)
                move_vcf(outputpath, outputpath + '.gz', 'text')
                os.remove(outputpath)
        # catches cases without genomic entries
        elif not self.hgvs_models or not self.get_variants():
            LOGGER.debug('VCF generation for case %s not possible, Error message: No variants',self.case_id)
        else:
            LOGGER.debug("Generating VCF for case %s", self.case_id)
            self.vcf = self.create_vcf(path)
            self.dump_vcf(path)
