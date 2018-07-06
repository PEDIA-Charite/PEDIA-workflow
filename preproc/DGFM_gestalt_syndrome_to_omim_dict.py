#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################################
# Parse 21147 from DGFM and output gestalt support omim
##########################################################

import json
import os
import csv 
import getopt
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-c', '--case', help='path to convert file')
    parser.add_argument('-o', '--output', help='path to output file')

    args = parser.parse_args()
    # syn_dict = {'name':[omim_id]}
    syn_dict = {}
    case_path = args.case
    
    #for case in os.listdir(case_path):
    file_name = '21147.json' 
    case_file = open(os.path.join(case_path, file_name))
    case_data = json.load(case_file)
    data = case_data['detected_syndromes_by_gestalt']
    for rank in data:
        syn_data = data[rank]
        name = syn_data['syndrome_name']
        omim = syn_data['omim_id']
        if type(name) == list:
            continue
        if syn_data['has_mask'] == 0:
            continue
        if name in syn_dict:
            continue
        omim_out = []
        if type(omim) == int:
            omim_out.append(omim)
        elif omim == None:
            omim_out = []
        else:
            omim_out = [omim_id for omim_id in omim if omim_id != -1]
        syn_dict.update({name:omim_out})
    print(len(syn_dict))
    print(syn_dict)

    syn_file = open(os.path.join(args.output,'gestalt_syn_to_omim.csv'), 'w')
    syn_writer = csv.writer(syn_file, delimiter='\t')
    for syn in syn_dict:
        syn_writer.writerow([syn, syn_dict[syn]])
    syn_file.close()
