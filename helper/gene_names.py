import os
import json
from lib.model.json import NewJson
from lib.model.case import Case
from lib.model.config import ConfigManager

config = ConfigManager()
# Initial Quality check of new json
unprocessed_jsons = os.path.join(config.aws['download_location'], 'cases')
case_id_list = [
    # INSERT CASE IDS HERE
]

still_empty = [
    # INSERT CASE IDS FOR SECOND PROCESSING
]



case_id_list = still_empty

json_files = [os.path.join(unprocessed_jsons, x)
              for x in os.listdir(unprocessed_jsons)
              if os.path.splitext(x)[1] == '.json'
              and int(os.path.splitext(x)[0]) in case_id_list]

corrected = config.preprocess['corrected_location']
new_json_objs = [NewJson.from_file(f, corrected) for f in json_files]
raw_dict = {}
for j in new_json_objs:
    print(j.get_case_id())
    print(j.get_js()['genomic_entries'])

# new_cases = [Case(j) for j in new_json_objs if j.check()[0]]

# create a list of case_id to gene_name
# case_json = [{c.case_id: j.entry_id} for c in new_cases for j in c.hgvs_models]
# print(case_json)
# case_json = [{c.case_id: j.get_json()} for c in new_cases for j in c.hgvs_models]
# json.dump(case_json, open('case_json_entries.json', 'w'), indent=4)
# gene_dict = {}
# case_ids = []
# case_ids_added = []
# for c in new_cases:
#     for j in c.hgvs_models:
#         case_ids.append(c.case_id)
#         if j.gene['gene_symbol']:
#             gene_dict[c.case_id] = j.gene['gene_symbol']
#             case_ids_added.append(c.case_id)
# case_ids = set(case_ids)
# case_ids_added = set(case_ids_added)
# case_diff = case_ids - case_ids_added
# print(case_diff)
# json.dump(gene_dict, open('gene_dict.json', 'w'), indent=4)
# print(gene_dict)

# case_gene_name = [{c.case_id: j.gene['gene_symbol']}
#                   for c in new_cases for j in c.hgvs_models]
# 
# empty_cases = [ c for c in case_gene_name if list(c.values())[0] == ""]
# print(empty_cases)
# # json.dump(case_gene_name, open('case_gene_name.json', 'w'))
