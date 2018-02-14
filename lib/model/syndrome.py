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

class SyndromeList:
    '''Structured object containing a number of syndromes for easier interaction.
    Data is saved in a dict by concatenating all omim keys.
    Secondary dicts can map additional information, such as single omim ids to these keys.

    Mappings between omim ids and syndromes is n:n.
    '''
    def __init__(self, syndrome_list=[]):
        self._data = {}
        if syndrome_list:
            self._save_list(syndrome_list)
            self._update_indexes()

    def add_syndrome(self, syndrome):
        for o in syndrome.omim_id:
            if o in self._mimdict:
                return False
        self._data[';'.join(syndrome.omim_id)] = syndrome
        self._update_indexes()
        return True

    def to_genes(self, omim):
        '''Convert syndromes to associated genes.
        '''
        pass

    def set_or_add_diagnosis(self, diagnosis):
        '''Add a diagnosis dictionary or flip confirmed on an existing diagnosis.
        '''
        r = self.get_syndrome_by_omim(diagnosis['omim_id'])
        if r:
            r.confirmed = True
        else:
            self.add_syndrome(Syndrome(**diagnosis))

    def get_diagnoses(self):
        return [ d for d in self._data.values() if d.confirmed ]

    def get_syndrome_by_omim(self,omim_id):
        if omim_id in self._mimdict:
            return self._data[self._mimdict[omim_id]]

    def _save_list(self, listdata):
        self._data = {";".join(l.omim_id) : l for l in listdata }

    def _update_indexes(self):
        '''Omim phenotypic ids are n:n mapped to syndromes.
        '''
        mimdict = {}
        for k,v in self._data.items():
            for o in v.omim_id:
                if o in mimdict:
                    mimdict[o].append(k)
                mimdict[o] = [k]
        self._mimdict = mimdict

    # overriding some python defaults for cleaner code structure
    def __contains__(self, omim_id):
        return bool(self.get_syndrome_by_omim(omim_id))

    def __getitem__(self, omim_id):
        return self.get_syndrome_by_omim(omim_id)

class Syndrome:
    '''
    Args:
        syndrome_name - Name of syndrome as string
        omim_id - Syndrome omim id
        has_mask - Whether the disease has a gestalt mask in Face2Gene
        combined_score - Feature + Gestalt score from F2G
        feature_score - Calculated from given features in F2G
        gestalt_score - from picture of case
        pheno_score - from phenomizer
        boqa_score - from boqa, another algorithm, based on the phenomizer
        confirmed - Whether the syndrome has been diagnosed
        gene_info - dict of gene-id to gene-symbol mapping
    '''
    def __init__(self, syndrome_name, omim_id, has_mask=False, combined_score=-1.0, feature_score=-1.0, gestalt_score=-1.0, pheno_score=-1.0, boqa_score=-1.0, confirmed=False, gene_info=[]):
        self.name = syndrome_name
        self.has_mask = has_mask
        # always coerce omim_id into a list
        omim_id = not isinstance(omim_id, list) and [omim_id] or omim_id
        self.omim_id = [ str(o) for o in omim_id ]
        self.combined_score = combined_score
        self.feature_score = feature_score
        self.gestalt_score = gestalt_score
        self.pheno_score = pheno_score
        self.boqa_score = boqa_score

        self.confirmed = confirmed

        self._genes = self._gene_info_to_genes(gene_info)

    def _gene_info_to_genes(self, gene_info):
        self._genes = gene_info

    def pheno_to_gene_dict(self, omim):
        genes = {}
        # load from gene info
        if self._gene_info:
            gene_mims = list(self._gene_info.keys())
            for gene_omim_id, gene_symbol in self._gene_info.items():
                gene_id = omim.mim_gene_to_entrez_id(gene_omim_id)
                genes[gene_omim_id] = {
                        'gene_id' : gene_id, 'gene_symbol':gene_symbol, 'gene_omim_id' : gene_omim_id
                        }
        # try to find additional information by searching phenotype to genotype mapping
        for o in self.omim_id:
            r = omim.mim_pheno_to_gene(o)
            # add genes not already mapped
            for k,v in r.items():
                if k not in genes:
                    genes[k] = v
        return genes

    def add_pheno_scores(self, pheno, boqa, gene_info=[]):
        '''Add phenomization and boqa scores.
        '''
        self.pheno_score = pheno
        self.boqa_score = boqa
        self._genes = gene_info or self._genes

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
        self.syndromes[syndrome.omim_id] = syndrome
        self.gene_omim_id = gene_omim_id
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol

    def add_syndrome(self, syndrome):
        if syndrome.omim_id not in self.syndromes:
            self.syndromes[syndrome.omim_id] = syndrome
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
