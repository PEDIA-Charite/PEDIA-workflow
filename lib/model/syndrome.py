'''
Syndrome model
---
Implements a model containing syndrome information.
This includes:
    omim_id - Identifier for Online Mendelian Inheritance in Man, is a list
    has_mask - Whether a mask in Face2Gene is available
    name - name of Syndrome
optionally:
    combined_score
    feature_score
    gestalt_score

    confirmed - whether the diagnosis has been confirmed, clinically or molecularly
        this information is currently approximated through the selected syndromes, interface
'''

def syndromes_to_genes(syndromes, omim):
    genes = {}
    for syndrome in syndromes:
        gene_list = syndrome.pheno_to_gene(omim)
        for g in gene_list:
            mim_id = g['gene_omim_id']
            if mim_id in genes:
                genes[mim_id].add_syndrome(syndrome)
            else:
                genes[mim_id] = Gene(**g, syndrome=syndrome)
    return genes

def update_genes(genes, pheno_boqa, omim):
    for pb in pheno_boqa.to_dict('records'):
        mim_id = omim.entrez_id_to_mim_gene(pb['gene-id'])
        if mim_id not in genes:
            s = Syndrome()
            genes[mim_id] = Gene(
                    gene_omim_id=mim_id,
                    gene_id=pb['gene-id'],
                    gene_symbol=pb['gene-symbol'],
                    syndrome=s)
        genes[mim_id].add_pheno_scores(pheno=pheno_boqa['value_pheno'],boqa=pheno_boqa['value_boqa'])
    return genes

class Syndrome:
    def __init__(self, syndrome_name, has_mask, omim_id, combined_score=-1.0, feature_score=-1.0, gestalt_score=-1.0, confirmed=False):
        self.name = syndrome_name
        self.has_mask = has_mask
        # always coerce omim_id into a list
        self.omim_id = not isinstance(omim_id, list) and [omim_id] or omim_id
        self.combined_score = combined_score
        self.feature_score = feature_score
        self.gestalt_score = gestalt_score

        self.confirmed = confirmed

    def pheno_to_gene(self, omim):
        res = []
        for o in self.omim_id:
            res += omim.mim_pheno_to_gene(o)
        return res

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __hash__(self):
        return hash(self.name.encode('utf-8'))

    def __repr__(self):
        return "{name} OMIM:{mim} Scores: feature:{fs} gestalt:{gs} combined:{cs}".format(
                name=self.name,
                mim=self.omim_id,
                fs=self.feature_score,
                gs=self.gestalt_score,
                cs=self.combined_score)

class Gene:
    '''Single mutated gene associated with a syndrome.
    For the candidate list we will need to infer this information from the omim_id to gene results.

    One phenotype might have multiple genetic associations. Try to first capture all of them here.
    '''
    def __init__(self, gene_omim_id, gene_id, gene_symbol, syndrome):
        if syndrome is None:
            raise RuntimeError(repr(syndrome))
        self._syndromes = {syndrome}
        self.gene_omim_id = gene_omim_id
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol

    def add_syndrome(self, syndrome):
        self._syndromes.add(syndrome)
        return self._syndromes

    def _pheno_to_gene(self, omim):
        res = []
        for o in self.omim_id:
            res += omim.mim_pheno_to_gene(o)
        return res

    def __hash__(self):
        return hash(self.gene_omim_id.encode('utf-8'))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
