import requests
import os
import argparse
import json
import csv


# Parse JSON content and convert
def convert_old(content, omim_dict):
    data = {}
    content = json.loads(content)
    case_data = content['case_data']
    data['case_id'] = case_data['case_id']
    # selected syndromes
    selected_syns = []
    for syn in case_data["selected_syndromes"]:
        selected_syn = {}
        selected_syn['diagnosis'] = syn['diagnosis']
        selected_syn['syndrome_name'] = syn['syndrome']['syndrome_name']
        selected_syn['omim_id'] = syn['syndrome']['omim_id']
        selected_syn['has_mask'] = 0
        selected_syns.append(selected_syn)
    data['selected_syndromes'] = selected_syns

    # features
    features = []
    for feature in case_data['selected_features']:
        features.append(feature['feature']['hpo_full_id'])
    data['features'] = features

    # detected syndromes
    detected_syns = []
    for syn in case_data['suggested_syndromes']:
        if 'syndrome' not in syn:
            continue
        detected_syn = {}
        detected_syn['syndrome_name'] = syn['syndrome']['syndrome_name']
        if syn['syndrome']['omim_id'] == None:
            continue
        if detected_syn['syndrome_name'] in omim_dict:
            detected_syn['omim_id'] = omim_dict[detected_syn['syndrome_name']]
        else:
            if len(syn['syndrome']['omim_ids']) == 0:
                if syn['syndrome']['omim_id'] == 0 or syn['syndrome']['omim_id'] == '0':
                    continue
                detected_syn['omim_id'] = syn['syndrome']['omim_id']
            else:
                detected_syn['omim_id'] = [omim for omim in syn['syndrome']['omim_ids'] if omim != 0]

        detected_syn['combined_score'] = 0
        detected_syn['gestalt_score'] = float(syn['gestalt_score'])
        detected_syn['feature_score'] = float(syn['feature_score'])
        detected_syn['has_mask'] = syn['syndrome']['app_valid']
        detected_syns.append(detected_syn)

    data['detected_syndromes'] = detected_syns
    data['algo_deploy_version'] = case_data['algo_version']
    data['genomic_entries'] = []
    sub = {}
    sub['user_email'] = "demo@pedia-study.org"#case_data['posting_user']['userEmail']
    sub['user_name'] = "pedia@study"#case_data['posting_user']['userDisplayName']
    sub['user_team'] = "pedia"#case_data['posting_user']['userInstitution']

    data['submitter'] = sub
    data['documents'] = [] if 'documents' not in content else content['documents']
    return data

# Parse input arguments
parser = argparse.ArgumentParser(description='Get case json file from lab')

parser.add_argument('-c', '--case', help='Select case id')
parser.add_argument('-o', '--output', help='Output folder')
parser.add_argument('-m', '--mapping', help='Omim mapping file 216_gestalt_syn_to_omim_final.csv')

args = parser.parse_args()
case_id = args.case
out_dir = args.output

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

omim_filename = args.mapping
omim_file = open(omim_filename, 'r')
reader = csv.reader(omim_file, delimiter='\t')
omim_dict = {}
for row in reader:
    if row[1] == "[]":
        continue
    omim = row[1][1:-1].split(', ')
    omim_dict[row[0]] = omim

case_file = open(case_id, 'r')
case_content = json.load(case_file)
output = convert_old(json.dumps(case_content), omim_dict)
out_file = open(os.path.join(out_dir, str(output['case_id']) + ".json"), 'w')
json.dump(output, out_file)
out_file.close()


