#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import csv 
import getopt
import sys
import argparse

######################################################################
# Parse all JSON files to construct syndrome to omim ID mapping        
#######################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-c', '--case', help='path to convert file')
    parser.add_argument('-o', '--output', help='path to output file')

    args = parser.parse_args()
    # syn_dict = {'name':[omim_id]}
    syn_dict = {}
    case_path = args.case
    for case in os.listdir(case_path):
        print(case)
        file_name = os.path.join(case_path, case)
        case_file = open(file_name)
        case_data = json.load(case_file)
        if 'detected_syndromes' not in case_data:
            continue
        data = case_data['detected_syndromes']
        for syn_data in data:
            name = syn_data['syndrome_name']
            omim = syn_data['omim_id']
            if type(name) == list:
                continue
            if syn_data['has_mask'] == 0:
                continue
            if name in syn_dict and omim == None:
                continue
            omim_out = []
            if type(omim) == int:
                omim_out.append(omim)
            else:
                omim_out = omim
            syn_dict.update({name:omim_out})
    print(len(syn_dict))
    print(syn_dict)

    syn_file = open(os.path.join(args.output,'gestalt_syn_to_omim.csv'), 'w')
    syn_writer = csv.writer(syn_file, delimiter='\t')
    for syn in syn_dict:
        print(syn)
        print(syn_dict[syn])
        syn_writer.writerow([syn, syn_dict[syn]])
    syn_file.close()
