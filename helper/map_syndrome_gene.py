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
    parser.add_argument('-o', '--output', help='path to output file')

    args = parser.parse_args()
    syn_path = args.syndrome
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
            gene = row[1].split(', ')[0]
            if omim not in omim_dict:
                omim_dict.update({omim:[gene]})
            else:
                if gene not in omim_dict[omim]:
                    omim_dict[omim].append(gene)
    gene_out = []
    print(omim_dict)
    omim_file.close()

    ############## Start mapping syndrome to gene
    output_filename = os.path.join(args.output,'gestalt_gene_phenotype_table.csv')
    out_file = open(output_filename, 'w')
    out_writer = csv.writer(out_file, delimiter='\t')

    #filename = '216_syndrome_omim.csv'
    filename = syn_path
    input_file = open(filename, 'r')
    reader = csv.reader(input_file, delimiter=',')
    for row in reader:
        tmp_gene = []
        omim = row[1]
        if '[' in omim:
            omim = omim.replace('[','')
            omim = omim.replace(']','')
            if omim == '' or omim == 0:
                continue
            omim = omim.split(', ')
            omim = [int(omim_id) for omim_id in omim]
        if omim == None or omim == 0:
            print(row[0])
        elif type(omim) == list:
            for omim_id in omim:
                if str(omim_id) in omim_dict:
                    for gene in omim_dict[str(omim_id)]:
                        if gene not in tmp_gene:
                            tmp_gene.append(gene)
                        if gene not in gene_out:
                            gene_out.append(gene)
        else:
            if str(omim) in omim_dict:
                for gene in omim_dict[str(omim)]:
                    if gene not in tmp_gene:
                        tmp_gene.append(gene)
                    if gene not in gene_out:
                        gene_out.append(gene)
            else:
                print(omim)
        out_writer.writerow([row[0], row[1], tmp_gene])
    out_file.close()

    print(len(gene_out))
    output_filename = os.path.join(args.output,'gestalt_gene_table.csv')
    out_file = open(output_filename, 'w')
    out_writer = csv.writer(out_file, delimiter='\t')
    for gene in gene_out:
        out_writer.writerow([gene])
    out_file.close()

