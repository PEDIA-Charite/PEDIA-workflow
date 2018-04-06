'''
Case model created from json files.
'''
import logging
from typing import Union, Dict

import pandas

from lib.model.json import OldJson, NewJson
from lib.utils import explode_df_column
from lib import constants


LOGGER = logging.getLogger(__name__)


def genes_to_single_cols(rowdata: pandas.Series) -> pandas.Series:
    '''Separate gene dict in genes column into separate columns.
    '''
    gene_dict = rowdata['genes']
    rowdata.drop('genes', inplace=True)
    genes = pandas.Series(gene_dict)
    rowdata = rowdata.append(genes)
    return rowdata


def create_gene_table(rowdata: pandas.Series, omim: 'Omim') -> pandas.Series:
    '''Get gene information from row information.
    This includes: all scores, gene_symbol, gene_id, gene_omim_id, syndrome_id
    '''
    disease_id = rowdata['omim_id']
    phenotypic_series = rowdata["phenotypic_series"]

    syndrome_name = rowdata["syndrome_name"] \
        or rowdata["disease-name_pheno"] \
        or rowdata["disease-name_boqa"]

    # get dict containing gene_id, gene_symbol, gene_omim_id
    genes = list(omim.mim_pheno_to_gene(disease_id).values())

    if rowdata["gene-id"]:
        genes += [
            {
                "gene_id": eid,
                "gene_symbol": omim.entrez_id_to_symbol(eid),
                "gene_omim_id": omim.entrez_id_to_mim_gene(eid)
            }
            for eid in rowdata["gene-id"].split(", ")
            if eid not in [g["gene_id"] for g in genes]
        ]

    # get all three scores provided by face2gene
    gestalt_score = rowdata["gestalt_score"]
    feature_score = rowdata['feature_score']
    combined_score = rowdata['combined_score']
    pheno_score = rowdata['value_pheno']
    boqa_score = rowdata['value_boqa']

    resp = pandas.Series({
        "disease_id": disease_id,
        "phenotypic_series": phenotypic_series,
        "syndrome_name": syndrome_name,
        "genes": genes,
        "gestalt_score": gestalt_score,
        "feature_score": feature_score,
        "combined_score": combined_score,
        "pheno_score": pheno_score,
        "boqa_score": boqa_score})
    return resp


def filter_phenotypic_series(ps_group: "DataFrame") -> "DataFrame":
    '''Filter phenotypic series group to return maximum values in group.'''
    if ps_group["phenotypic_series"].iloc[0] == "":
        return ps_group

    ps_group = ps_group.groupby("gene_id", as_index=False).max()
    return ps_group


class Case:
    '''
    Exposes the following properties:
    case_id - Unique identifier for the case in Face2Gene
    variants - list of hgvs objects describing valid hgvs variants
    analyzed_syndromes - dictionary of syndromes using omim id as key
    features - list of hpo terms
    diagnosis - list of syndromes selected as diagnosis
    submitter - submitter information containing fields for email, name, team
    vcf - list of vcf filenames
    '''

    def __init__(self, data: Union[OldJson, NewJson],
                 error_fixer: "ErrorFixer",
                 exclude_benign_variants: bool = True):
        self.algo_version = data.get_algo_version()
        self.case_id = data.get_case_id()

        # get both the list of hgvs variants and the hgvs models used in the
        # parsing
        self.hgvs_models = data.get_variants(error_fixer)
        self.syndromes = data.get_syndrome_suggestions_and_diagnosis()
        self.features = data.get_features()
        self.submitter = data.get_submitter()
        self.vcf = data.get_vcf()
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

        # fill nans created by merge
        self.syndromes.fillna(
            {
                'combined_score': 0.0,
                'feature_score': 0.0,
                'gestalt_score': 0.0,
                'value_pheno': 0.0,
                'value_boqa': 0.0,
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

    def get_gene_list(self, omim: 'Omim', recreate: bool = False,
                      filter_entrez_id: bool = True) \
            -> Dict[str, str]:
        '''Get list of genes from syndromes. Save them back to self.genes
        for faster lookup afterwards.
        Args:
            filter_entrez_id: Only include genes with existing entrez gene ids.
        '''
        # return existing gene list, if it already exists
        if self.gene_scores is not None and not recreate:
            return self.gene_scores

        LOGGER.debug("Generating geneList for case %s", self.case_id)

        # add or update phenotypic series information to syndromes table
        self.syndromes["phenotypic_series"] = \
            self.syndromes["omim_id"].astype(str).apply(
                omim.omim_id_to_phenotypic_series
            )

        # # group by phenotypic series and return highest unless empty
        # syndrome_series = self.syndromes.groupby("phenotypic_series").apply(
        #     filter_phenotypic_series
        # )
        # syndrome_series = syndrome_series.reset_index(drop=True)

        gene_table = self.syndromes.apply(
            lambda x: create_gene_table(x, omim), axis=1
        )

        # explode the gene table on genes to separate the genetic entries
        gene_table = explode_df_column(gene_table, 'genes')
        gene_table = gene_table.apply(genes_to_single_cols, axis=1)

        # only select entries with non-empty gene ids
        if filter_entrez_id:
            gene_table = gene_table.loc[gene_table["gene_id"] != ""]

        # reset indexing for grouping
        gene_table = gene_table.reset_index(drop=True)

        # group by phenotypic series and return highest unless empty
        gene_table = gene_table.groupby("phenotypic_series").apply(
            filter_phenotypic_series
        )

        gene_scores = gene_table.to_dict('records')
        self.gene_scores = gene_scores
        return gene_scores

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
        # check maximum gestalt score
        max_gestalt_score = max(self.syndromes['gestalt_score'])
        if max_gestalt_score <= 0:
            issues.append(
                {
                    "type": "NO_GESTALT",
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
                omim.omim_id_to_phenotypic_series(d)
                for d in diagnosis["omim_id"]
            ]
            if len(set(diagnosis_series)) != 1:
                issues.append(
                    {
                        "type": "MULTI_DIAGNOSIS",
                        "data": list(diagnosis["omim_id"])
                    }
                )
                valid = False
        else:
            LOGGER.warning("No omim object. Some checks will not run.")

        if not self.get_variants():
            issues.append(
                {
                    "type": "NO_GENOMIC"
                }
            )
            valid = False

        # check for benign exclusion
        issues.append(
            {
                "type": "REMOVE_NORMAL_VARIANTS",
                "data": {
                    "all": len(self.get_variants(exclusion=False)),
                    "excluded": len(self.get_variants(exclusion=True))
                }
            }
        )

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
        return not self.vcf
