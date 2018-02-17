import csv
import json

FILE_PATH = '/data/users/Max/PEDIA-workflow/handcurated_cases.csv'
JSON_PATH = '/data/users/Max/PEDIA-workflow/hgvs_errors.json'

corrected_entries = {}
with open(FILE_PATH, 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file, delimiter=',', quotechar='"')
    prev_id = None
    for row in csv_reader:
        mutation_id = row['genomic entry']
        # multiple entries are written without mutation_id
        if not mutation_id:
            mutation_id = prev_id
        prev_id = mutation_id
        # remove trailing whitespace
        corrected_hgvs_code = row['correct HGVS'].strip()
        if not corrected_hgvs_code:
            continue


        if mutation_id in corrected_entries:
            corrected_entries[mutation_id].append(corrected_hgvs_code)
        else:
            corrected_entries[mutation_id] = [corrected_hgvs_code]

hgvs_errordict = json.load(open(JSON_PATH, 'r'))

for entry_id, corrected_hgvs in corrected_entries.items():
    print(entry_id)
    print(corrected_hgvs)
    hgvs_errordict[entry_id]['cleaned'] = corrected_hgvs

json.dump(hgvs_errordict, open('hgvs_errors_load_handcrafted.json', 'w'),
          indent=4)
