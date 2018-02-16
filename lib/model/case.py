'''
Case model created from json files.
'''
from functools import reduce
from typing import Union, Dict

import pandas

from lib.model.json import OldJson, NewJson


def create_gene_table(rowdata: pandas.Series, omim: 'Omim') -> pandas.Series:
    '''Get gene information from row information.
    This includes: all scores, gene_symbol, gene_id, gene_omim_id, syndrome_id
    '''
    disease_id = rowdata['omim_id']
    # get dict containing gene_id, gene_symbol, gene_omim_id
    genes = list(omim.mim_pheno_to_gene(disease_id).values())
    # get all three scores provided by face2gene
    gestalt_score = pandas.notna(rowdata['gestalt_score']) \
        and rowdata['gestalt_score'] or 0.0
    feature_score = pandas.notna(rowdata['gestalt_score']) \
        and rowdata['feature_score'] or 0.0
    combined_score = pandas.notna(rowdata['combined_score']) \
        and rowdata['combined_score'] or 0.0
    # get scores from phenomizer
    if 'value_pheno' in rowdata:
        pheno_score = pandas.notna(rowdata['value_pheno']) \
            and rowdata['value_pheno'] or 0.0
    else:
        pheno_score = 0.0

    if 'value_boqa' in rowdata:
        boqa_score = pandas.notna(rowdata['value_boqa']) \
            and rowdata['value_boqa'] or 0.0
    else:
        boqa_score = 0.0

    resp = pandas.Series({
        "disease_id": disease_id,
        "genes": genes,
        "gestalt_score": gestalt_score,
        "feature_score": feature_score,
        "combined_score": combined_score,
        "pheno_score": pheno_score,
        "boqa_score": boqa_score})
    return resp


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

    def __init__(self, json_object: Union[OldJson, NewJson]):
        self._from_json_object(json_object)

    def _from_json_object(self, data: Union[OldJson, NewJson]):
        '''Load case information from json object.
        '''
        self.case_id = data.get_case_id()
        self.variants = data.get_variants()
        self.syndromes = data.get_syndrome_suggestions_and_diagnosis()
        self.features = data.get_features()
        self.submitter = data.get_submitter()
        self.vcf = data.get_vcf()

    def phenomize(self, pheno: 'PhenomizerService') -> bool:
        '''Add phenomization information to genes from boqa and phenomizer.
        Args:
            omim: Omim object to handle id translation.
            pheno: PhenomizerService to handle API calls for phenomizer and
                   boqa.
        '''
        pheno_boqa = pheno.disease_boqa_phenomize(self.features)
        if pheno_boqa is None:
            return False
        pheno_boqa.index = pheno_boqa.index.astype(int)
        # merge pheno and boqa scores dataframe with our current syndromes
        # dataframe which contains face2gene scores
        self.syndromes = self.syndromes.merge(
            pheno_boqa, left_on='omim_id', how='outer', right_index=True)
        return True

    def get_gene_list(self, omim: 'Omim', recreate: bool=False) \
            -> Dict[str, str]:
        '''Get list of genes from syndromes. Save them back to self.genes
        for faster lookup afterwards.
        '''
        # return existing gene list, if it already exists
        if hasattr(self, 'gene_scores') and not recreate:
            if self.gene_scores:
                return self.gene_scores

        gene_table = self.syndromes.apply(
            lambda x: create_gene_table(x, omim), axis=1)
        gene_scores = gene_table.to_dict('records')
        self.gene_scores = gene_scores
        return gene_scores

    def check(self) -> bool:
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
        max_gestalt_score = reduce(
            lambda x, y: max(x, y.gestalt_score), self.detected_syndromes, 0.0)
        if max_gestalt_score <= 0:
            issues.append('Maximum gestalt score is 0. \
                          Probably no image has been provided.')
            valid = False

        # check that only one syndrome has been selected
        if len(self.diagnosis) != 1:
            issues.append(
                '{} syndromes have been selected. Only 1 syndrome should be \
                selected for PEDIA inclusion.'.format(len(self.diagnosis)))
            valid = False

        # check that molecular information is available at all
        if len(self.variants) == 0:
            issues.append('No valid genomic entries available.')
            valid = False

        return valid, issues

    def eligible_training(self) -> bool:
        '''Eligibility of case for training.
        Exclusion criteria are:
            available real vcf
        '''
        return not self.vcf
