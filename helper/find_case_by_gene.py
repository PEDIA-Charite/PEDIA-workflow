#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import csv 
import getopt
import sys
import argparse

def get_case_info(file_writer, file_content, reason, status="Pass"):
    tmp_out = []
    info = []
   
    genomic = file_content['genomicData']
    ge_out = [ge_entry['Test Information']['Gene Name'] for ge_entry in genomic]
    for gene in gene_list:
        if gene in ge_out:
            # Output results to csv
            tmp_out.append(fileName.split('.')[0])
            tmp_out.append(gene)
            out.append(tmp_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-c', '--case', help='path to convert file')
    parser.add_argument('-o', '--output', help='path to output file')
    parser.add_argument('-l', '--gene-list', help='list of genes you want to find in case')
    
    args = parser.parse_args()
    case_path = args.case
    list_file = open(args.gene_list, 'r')
    reader = csv.reader(list_file)
    gene_list = []
    out = []
    for row in reader:
        gene_list.append(row[0])
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    output_filename = os.path.join(args.output, 'cases.csv')
    file_in_dir = [case for case in os.listdir(case_path) if case.split('.')[-1] == 'json']
    with open(output_filename, 'w') as csvfile:
        file_writer = csv.writer(csvfile, delimiter='\t')

        for fileName in file_in_dir:
            test = fileName
            print(test)
            file_content = json.load(open(os.path.join(case_path, fileName)))
            get_case_info(file_writer, file_content, [])

        s = sorted(out, key = lambda x: (x[1], int(x[0])))
        for row in s:
            file_writer.writerow(row)
