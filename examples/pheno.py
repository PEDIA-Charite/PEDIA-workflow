import os
import sys
import pickle
sys.path.append(os.getcwd())
from lib.model.config import ConfigManager
from lib.api.phenomizer import PhenomizerService
from lib.api.omim import Omim


cm = ConfigManager()

pheno = PhenomizerService(config=cm)
omim = Omim(config=cm)

# load cases with previous preprocessing applied for faster testing
case_list = pickle.load(open('case_cleaned.p', 'rb') )

t = case_list[0]

# r = pheno.disease_boqa_phenomize(t.features)
for c in case_list:
    c.phenomize(omim, pheno)
