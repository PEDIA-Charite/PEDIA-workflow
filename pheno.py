import pickle
from lib.model.config import ConfigManager
from lib.api.phenomizer import PhenomizerService
from lib.api.omim import Omim
from lib.model.syndrome import update_genes

cm = ConfigManager()

pheno = PhenomizerService(config=cm)
omim = Omim(config=cm)

case_list = pickle.load(open('case_cleaned.p', 'rb') )

t = case_list[0]

# pb = pheno.get_phenomize(t.features)
pb = pheno.pheno_boqa(t.features, omim)
print(pb)
# update_genes(t.genes, pb, omim)
