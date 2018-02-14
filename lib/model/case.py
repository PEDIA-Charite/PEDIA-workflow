from lib.model.json import OldJson, NewJson
from lib.model.syndrome import Gene, Syndrome

'''
Case model created from json files.
'''

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

    def __init__(self, json_object):
        if isinstance(json_object, OldJson):
            self._from_old_json(json_object)
        elif isinstance(json_object, NewJson):
            self._from_new_json(json_object)
        else:
            raise TypeError('Unsupported Type')

    def _from_old_json(self, data):
        raise RuntimeError('This has not been implemented')

    def _from_new_json(self, data):
        self.case_id = data.get_case_id()
        self.variants = data.get_variants()
        self.syndromes = data.get_syndrome_suggestions_and_diagnosis()
        self.features = data.get_features()
        self.submitter = data.get_submitter()
        self.vcf = data.get_vcf()

    def syndromes_to_genes(self, omim):
        return self.syndromes.to_genes(omim)

    def phenomize(self, pheno):
        '''Add phenomization information to genes from boqa and phenomizer.
        Args:
            omim - Omim object to handle id translation.
            pheno - PhenomizerService to handle API calls for phenomizer and boqa.
        '''
        pheno_boqa = pheno.disease_boqa_phenomize(self.features)
        if pheno_boqa is None:
            return False
        pheno_boqa.index = pheno_boqa.index.astype(int)
        r = self.syndromes.merge(pheno_boqa, left_on='omim_id', how='outer', right_index=True)
        self.syndromes = r
        return True

    def check(self):
        '''Check whether Case fulfills all provided criteria.

        The criteria are:
            picture has been provided - gestalt_score in detected_syndromes should be greater than 0
            clinical diagnosis - selected_syndromes should not be empty
            single monogenetic disease - not multiple syndromes selected and not multiple pathogenic mutations in different genes
            SNP mutations - no microdel/dup or other large scale aberrations

        Note that this function is more definitive than the json level check, as the validity of hgvs parsing has already been established.
        '''
        valid = True
        issues = []
        # check maximum gestalt score
        max_gestalt_score = functools.reduce(lambda x,y: max(x,y.gestalt_score), self.detected_syndromes, 0.0)
        if max_gestalt_score <= 0:
            issues.append('Maximum gestalt score is 0. Probably no image has been provided.')
            valid = False

        # check that only one syndrome has been selected
        if len(self.diagnosis) != 1:
            issues.append('{} syndromes have been selected. Only 1 syndrome should be selected for PEDIA inclusion.'.format(len(self.diagnosis)))
            valid = False

        # check that molecular information is available at all
        if len(self.variants) == 0:
            issues.append('No valid genomic entries available.')
            valid = False

        return valid, issues

    def eligible_training(self):
        '''Eligibility of case for training.
        Exclusion criteria are:
            available real vcf
        '''
        return not self.vcf
