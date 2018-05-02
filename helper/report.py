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

    # Check selected syndrome
    if 'selected_syndromes' in file_content:
        syn_out = []
        omim_out = []
        for syn in file_content['selected_syndromes']:
            syn_out.append(syn['syndrome_name'])
            omim_out.append(syn['omim_id'])

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
    tmp_out.append(omim_out)
    file_writer.writerow(tmp_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-c', '--case', help='path to convert file')
    parser.add_argument('-l', '--log', help='path to log file')
    parser.add_argument('-o', '--output', help='path to output file')
    parser.add_argument('-j', '--json-check', help='path to json_check_failed.log')
    
    args = parser.parse_args()
    case_path = args.case
    log_file = open(args.log)
    log_data = json.load(log_file)

    log_precheck_file = open(args.json_check)
    log_precheck_data = json.load(log_precheck_file)

    log_failed = [case for case in log_precheck_data if case["valid"] == False]
    output_filename = os.path.join(args.output, 'summary_cases.csv')
    
    # ['failed', 'benign_excluded', 'pathogenic_missing', 'vcf_failed', 'passed']
    with open(output_filename, 'w') as csvfile:
        file_writer = csv.writer(csvfile, delimiter=' ')

        for fileName in log_data["passed"]:
            test = fileName
            file_content = json.load(open(os.path.join(case_path, fileName + '.json')))
            get_case_info(file_writer, file_content, [])

        failed_data = log_data["failed"]
        for fileName in failed_data.keys():
            file_content = json.load(open(os.path.join(case_path, fileName + '.json')))
            reason = []
            for data in failed_data[fileName][1:][0]:
                reason.append(data['type'])
            get_case_info(file_writer, file_content, reason, "Fail")

        failed_data = log_data["vcf_failed"]
        for fileName in failed_data.keys():
            file_content = json.load(open(os.path.join(case_path, fileName + '.json')))
            reason = []
            for data in failed_data[fileName][1:][0]:
                reason.append(data['type'])
            get_case_info(file_writer, file_content, reason, "Fail")
            
        for case in log_failed:
            file_content = json.load(open(os.path.join('process/aws_dir/cases/', case['case_id'] + '.json')))
            fileName = case['case_id']
            reason = []
            for data in case['issues']:
                reason.append(data)
            get_case_info(file_writer, file_content, reason, "Fail")
