#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import csv 
import getopt
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-s', '--syndrome', help='path to convert file')
    parser.add_argument('-p', '--pedia', help='path to convert file')
    parser.add_argument('-o', '--output', help='path to output file')

    args = parser.parse_args()
    syn_path = args.syndrome
    pedia_path = args.pedia
    deprecated_file = open('data/omim_deprecated_replacement.json', 'r')
    deprecated_omim = json.load(deprecated_file)

    # Parse omim phenotype to gene mapping file
    omim_dict = {} 
    omim_file = open('data/morbidmap.txt', 'r')
    omim_reader = csv.reader(omim_file, delimiter='\t')
    count = 0
    for row in omim_reader:
        count = count + 1
        if len(row) > 2:
            omim = row[0].split(', ')[-1].split(' ')[0]
            gene = row[1].split(', ')[0][1:-1]
            if omim not in omim_dict:
                omim_dict.update({omim:[gene]})
            else:
                if gene not in omim_dict[omim]:
                    omim_dict[omim].append(gene)
    gene_out = []
    print(omim_dict)
    omim_file.close()

    ############## Start mapping syndrome to gene
    #output_filename = os.path.join(args.output,'pedia_gestalt_gene_phenotype_table.csv')
    #out_file = open(output_filename, 'w')
    #out_writer = csv.writer(out_file, delimiter='\t')

    #filename = '216_syndrome_omim.csv'
    filename = syn_path
    input_file = open(filename, 'r')
    reader = csv.reader(input_file, delimiter='\t')
    syns = []
    for row in reader:
        syns.append(row[0])
    input_file.close()

    filename = pedia_path
    input_file = open(filename, 'r')
    reader = csv.reader(input_file, delimiter='\t')
    p_syns = []
    p_no_syns = []
    print(syns)
    for row in reader:
        syn = row[7]
        syn = syn.replace('[','')
        syn = syn.replace(']','')
        syn = syn.split(', ')[0][1:-1]
        if syn in syns:
            if syn not in p_syns:
                p_syns.append(syn)
        else:
            if syn not in p_no_syns:
                p_no_syns.append(syn)
    print(len(p_syns))
    print(p_no_syns)
    input_file.close()

