import os
import csv 
import getopt
import sys
import argparse
import shutil
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-i', '--case', help='path to convert file')
    parser.add_argument('-f', '--f2g', help='path to f2g dgfm file')
    parser.add_argument('-s', '--syn-omim', help='path to f2g dgfm file')
    parser.add_argument('-o', '--output', help='path to output file')
    parser.add_argument('-l', '--list', help='path to config_gestalt.csv')

    # -s 216_gestalt_syn_to_omim_final.csv

    args = parser.parse_args()
    case_path = args.case
    dgfm_path = args.f2g
    omim_path = args.syn_omim
    name = args.list
    g_file = open(name, 'r')
    reader = csv.reader(g_file)
    file_list = []
    for row in reader:
        file_list.append(row[0] + '.json')
    g_file.close()

    omim_file = open(omim_path, 'r')
    reader = csv.reader(omim_file, delimiter='\t')
    omim = {}
    for row in reader:
        if row[1] == "[]":
            continue
        omim_id = row[1][1:-1].split(", ")
        omim[row[0]] = omim_id
    omim_file.close()

    for case in file_list:
        name = os.path.join(dgfm_path, case)
        content = json.load(open(name))
        syn_dict = {}

        # Parse Gestalt
        g_data = content['detected_syndromes_by_gestalt']
        for i in g_data:
            data = g_data[i]
            syn_id = data["syndrome_id"]
            name = data["syndrome_name"] 
            if name in omim:
                data['omim_id'] = omim[name]
            data['combined_score'] = 0
            syn_dict[syn_id] = data
        # Parse FM
        f_data = content['detected_syndromes_by_fm']
        syn_array = []
        for i in f_data:
            data = f_data[i]
            syn_id = data["syndrome_id"]
            if syn_id in syn_dict:
                continue
            name = data["syndrome_name"]
            if type(name) == list:
                name = str(name)
            data['syndrome_name'] = name
            data['combined_score'] = 0
            syn_dict[syn_id] = data
        for i in syn_dict:
            syn_array.append(syn_dict[i])

        pedia_name = os.path.join(case_path, case)
        pedia_content = json.load(open(pedia_name))
        pedia_content['detected_syndromes'] = syn_array

        if case == "249639.json" or case == "245514.json":
            selected_syn = content['selected_syndromes']
            selected_syn[0]['diagnosis'] = "MOLECULARLY_DIAGNOSED"
            selected_syn[0]['syndrome_name'] = selected_syn[0]['syndrome_name'][0]
            pedia_content['selected_syndromes'] = selected_syn
        if case == "226359.json" or case == "226371.json" or case == "226376.json" or case == "226358.json":
            selected_syn = pedia_content['selected_syndromes']
            selected_syn[0]['has_mask'] = 1
            selected_syn[0]['syndrome_name'] = "Hutchinson-Gilford Progeria Syndrome; HGPS"
            selected_syn[0]['omim_id'] = [176670]
            pedia_content['selected_syndromes'] = selected_syn

        if not os.path.exists(args.output):
            os.makedirs(args.output)
        out_name = os.path.join(args.output, case)
        with open(out_name, 'w') as outfile:
            json.dump(pedia_content, outfile)



