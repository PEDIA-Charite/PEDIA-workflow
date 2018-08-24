import requests
import os
import argparse
import json
import csv


# Parse JSON content and convert
def convert_old(content, omim_dict):
    data = {}
    content = json.loads(content)
    data['case_id'] = content['f2g_case_id']
    data['lab_case_id'] = content['lab_case_id']
    data['lab_id'] = content['lab_id']
    data['lab_id'] = content['lab_id']
    case_data = content['case_data']
    
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
    data['lab_info'] = content['lab_info']
    data['genomic_entries'] = []
    sub = {}
    sub['user_email'] = case_data['posting_user']['userEmail']
    sub['user_name'] = case_data['posting_user']['userDisplayName']
    sub['user_team'] = case_data['posting_user']['userInstitution']

    data['submitter'] = sub
    data['documents'] = []
    data['gender'] = case_data['gender']
    data['ethnicity'] = case_data['ethnicity_array']
    return data




# Parse input arguments

parser = argparse.ArgumentParser(description='Get case json file from lab')

parser.add_argument('-l', '--lab', help='Select lab number')
parser.add_argument('-c', '--case', help='Select case id')
parser.add_argument('-o', '--output', help='Output folder')
parser.add_argument('-m', '--mapping', help='Omim mapping file 216_gestalt_syn_to_omim_final.csv')

args = parser.parse_args()
lab = args.lab
case_id = args.case
out_dir = args.output

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

if lab == 'bonn':
    lab_id = '8'
    lab_name: 'bonn'
    api_key = "061fdb4a0efd8ede403d458e18a7cae2"
    api_secret = "Fi4B3PZWUH63wLLRvA9qLRzLKKMm8S"
elif lab == 'pedia':
    lab_id = 26
    lab_name =  "PEDIA"
    api_key = "380dfad0293b3803b49f1c7121736e46"
    api_secret = "2OfrqOCCUmsad0LEPI42Ip0BtiyYu8G"

lab_connect = "https://app.face2gene.com/api/labs/auth/connect"
lab_get_case_list = "https://app.face2gene.com/api/lab-cases?lab_id=" + lab_id
lab_get_case = "https://app.face2gene.com/api/lab-cases/" + case_id
response = requests.post(lab_connect, data=json.dumps({"api_key": api_key, "secret": api_secret}), headers={"Content-Type": "application/json", "Accept": "application/json"})

token = ''
if response.status_code == 200:
    token = json.loads(response.content.decode('utf-8'))['jwt']
else:
    print(response.status_code)

omim_filename = args.mapping
omim_file = open(omim_filename, 'r')
reader = csv.reader(omim_file, delimiter='\t')
omim_dict = {}
for row in reader:
    if row[1] == "[]":
        continue
    omim = row[1][1:-1].split(', ')
    omim_dict[row[0]] = omim 

if token != '':
    auth = "Bearer " + token
    response = requests.get(lab_get_case, headers={"Accept":"application/json", "Authorization":auth})
    if response.status_code == 200:
        case_content = json.loads(response.content.decode('utf-8'))
        output = convert_old(json.dumps(case_content), omim_dict)
        out_file = open(os.path.join(out_dir, str(output['case_id']) + ".json"), 'w')
        json.dump(output, out_file)
        out_file.close()
    else:
        print(response.status_code)


