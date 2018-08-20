#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import csv 
import getopt
import sys
import argparse

def get_case_info(file_writer, file_content, reason, status="Pass"):
    out = []
    tmp_out = []
    info = []
    ge_out = []
    hgvs_out = []
   
    # Submitter
    submitter = file_content["submitter"]
    name = submitter["user_name"]
    email = submitter["user_email"]
    team = submitter["user_team"]

    # Check disease-causing gene
    for entry in file_content['genomic_entries']:
        if type(entry) != int and 'variants' in entry:
            if 'notes' in entry['variants']:
                info.append(entry['variants']['notes'])
    if 'genomicData' in file_content:
        genomic = file_content['genomicData']
        ge_out = [ge_entry['Test Information']['Gene Name'] for ge_entry in genomic]
        hgvs_out = [ge_entry['Mutations']['HGVS-code'] for ge_entry in genomic]

    # Check selected syndrome
    if 'selected_syndromes' in file_content:
        syn_out = []
        omim_out = []
        diag_out = []
        for syn in file_content['selected_syndromes']:
            syn_out.append(syn['syndrome_name'])
            omim_out.append(syn['omim_id'])
            diag_out.append(syn['diagnosis'])

        if len(syn_out) > 1:
            list_syn = []
            single_syn = ""
            for syn in omim_out:
                if type(syn) == list:
                    list_syn = syn
                else:
                    single_syn = syn

            if single_syn not in list_syn:
                if status == "Pass":
                    print(fileName + ": Multiple different selected syndromes")

        elif len(syn_out) == 0:
            if status == "Pass":
                print(fileName + ' has no syndrome')

    # Output results to csv
    tmp_out.append(fileName.split('.')[0])
    tmp_out.append(status)
    tmp_out.append(str(reason)[1:-1])
    tmp_out.append(str(ge_out)[1:-1])
    tmp_out.append(name)
    tmp_out.append(email)
    tmp_out.append(team)
    tmp_out.append(syn_out)
    tmp_out.append(diag_out)
    tmp_out.append(omim_out)
    tmp_out.append(hgvs_out)
    file_writer.writerow(tmp_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-c', '--case', help='path to convert file')
    parser.add_argument('-l', '--log', help='path to log file')
    parser.add_argument('-o', '--output', help='path to output file')
    
    args = parser.parse_args()
    case_path = args.case
    log_file = open(args.log)
    log_data = json.load(log_file)

    name = 'config_gestalt.csv'
    g_file = open(name, 'r')
    reader = csv.reader(g_file)
    file_list = []
    for row in reader:
        file_list.append(int(row[0]))

    output_filename = os.path.join(args.output, 'summary_gestalt_cases.csv')
    
    # ['failed', 'benign_excluded', 'pathogenic_missing', 'vcf_failed', 'passed']
    with open(output_filename, 'w') as csvfile:
        file_writer = csv.writer(csvfile, delimiter='\t')

        for fileName in log_data["passed"].keys():
            if int(fileName) not in file_list:
                continue
            test = fileName
            file_content = json.load(open(os.path.join(case_path, fileName + '.json')))
            get_case_info(file_writer, file_content, [])

