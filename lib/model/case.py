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
from lib import vcf_jannovar
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

        self.vcf = None

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
                'differential': False,
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
            if "gene-id" in syndrome and syndrome["gene-id"]:
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

                # uniqueness constraint on disease_id and
                # gene_id
                key = "{}|{}".format(syndrome_name, gene["gene_id"])
                pheno_score = syndrome["pheno_score"] \
                    if "pheno_score" in syndrome else 0.0
                boqa_score = syndrome["boqa_score"] \
                    if "boqa_score" in syndrome else 0.0
                update_data = dict(
                    {
                        "disease_id": disease_id,
                        "phenotypic_series": phenotypic_series,
                        "syndrome_name": syndrome_name,
                        "gestalt_score": syndrome["gestalt_score"],
                        "feature_score": syndrome["feature_score"],
                        "combined_score": syndrome["combined_score"],
                        "pheno_score": pheno_score,
                        "boqa_score": boqa_score
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

        # check if there is at least one feature (HPO)
        features = self.features
        if len(features) < 1:
            issues.append(
                {
                    "type": "NO_FEATURES",
                }
            )
            valid = False

        scores = [
            "gestalt_score"
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
            ps_dict = {}
            for _, diag in diagnosis.iterrows():
                ps_res = omim.omim_id_to_phenotypic_series(
                    str(diag["omim_id"])
                ) or str(diag["omim_id"])
                # ignore entries without omim id
                if ps_res == '0':
                    continue
                if diag["syndrome_name"] in ps_dict:
                    ps_dict[diag["syndrome_name"]].add(ps_res)
                else:
                    ps_dict[diag["syndrome_name"]] = set([ps_res])
            # compact ps_dict based on omim ids
            reduced_ps_dict = {}

            for key, series in ps_dict.items():
                contained = False
                for other_key, other_series in ps_dict.items():
                    if other_key == key:
                        continue
                    if series <= other_series:
                        contained = True
                if not contained:
                    reduced_ps_dict[key] = series

            if len(reduced_ps_dict) > 1:
                issues.append(
                    {
                        "type": "MULTI_DIAGNOSIS",
                        "data": {
                            "orig": list(diagnosis["omim_id"]),
                            'names': list(reduced_ps_dict.keys()),
                            "converted_ids": [
                                e for v in reduced_ps_dict.values() for e in v
                            ]
                        }
                    }
                )
                valid = False
        else:
            LOGGER.warning("No omim object. Some checks will not run.")

        if not self.get_variants():
            raw_entries = self.data.get_genomic_entries()
            if not raw_entries:
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
        else:
            # Check if there are multiple different disease-causing genes 
            raw_entries = self.data.get_genomic_entries()
            entries = self.get_hgvs_models()
            if len(entries) > 1:
                genes = [entry.gene['gene_id'] for entry in entries if entry.gene['gene_id']]
                if len(set(genes)) > 1:
                    issues.append(
                        {
                            "type": "MULTI_DIFFERENT_DISEASE_CAUSING_GENE",
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

    def create_vcf_from_hgvs(
            self,
            hgvs_strings: [str],
            vcf_path: str,
            tmp_path: str,
    ) -> pandas.DataFrame:
        '''Generates vcf dataframe. If an error occurs the error message is
        returned.'''
        vcf_data = vcf_jannovar.create_vcf(
            hgvs_strings,
            self.get_hgvs_models()[0].zygosity.lower(),
            self.case_id,
            tmp_path
        )

        if isinstance(vcf_data, pandas.DataFrame):
            vcf_jannovar.write_vcfdf(vcf_data, vcf_path)
        elif isinstance(vcf_data, str):
            LOGGER.debug(
                "%s: VCF generation failed. Error: %s",
                self.case_id, vcf_data
            )
        else:
            LOGGER.debug(
                "%s: No vcf generated yet.",
                self.case_id
            )
        return vcf_data

    def put_hgvs_vcf(
            self,
            path: str,
            recreate: bool = False,
    ) -> None:
        '''Dumps vcf file to given path as <case_id>.vcf.gz.'''
        if not self.get_variants():
            LOGGER.debug(
                '%s: VCF generation impossible. Error: No variants',
                self.case_id
            )
            return

        hgvs_strings = [str(v) for v in self.get_variants()]
        vcf_path = os.path.join(path, self.case_id + ".vcf.gz")

        if not recreate and os.path.exists(vcf_path):
            vcf_data = vcf_jannovar.read_vcfdf(vcf_path)
            vcf_hgvs = vcf_jannovar.get_hgvs_codes(vcf_data)

            # ensure that hgvs strings are same
            if all(h in vcf_hgvs for h in hgvs_strings) \
                    and all(h in hgvs_strings for h in vcf_hgvs):
                LOGGER.debug("%s: Use existing vcf.", self.case_id)
        else:
            vcf_data = self.create_vcf_from_hgvs(
                hgvs_strings, vcf_path, path
            )
        self.vcf = vcf_data
