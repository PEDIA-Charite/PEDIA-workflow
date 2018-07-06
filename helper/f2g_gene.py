#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#######################################################################
# Return number of gene mapped from f2g gestalt syndrome
#######################################################################


import json
import os
import csv 
import getopt
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-o', '--output', help='path to output file')

    log_file = open("21147.json")
    case_data = json.load(log_file)

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
    g_syndrome = case_data['detected_syndromes_by_gestalt']
    for syn in g_syndrome:
        omim = g_syndrome[syn]['omim_id']
        if omim == None or omim == 0:
            print(g_syndrome[syn]['syndrome_name'])
        elif type(omim) == list:
            for omim_id in omim:
                if str(omim_id) in omim_dict:
                    for gene in omim_dict[str(omim_id)]:
                        if gene not in gene_out:
                            gene_out.append(gene)
        else:
            if str(omim) in omim_dict:
                for gene in omim_dict[str(omim)]:
                    if gene not in gene_out:
                        gene_out.append(gene)
            else:
                print(omim)

    print(len(gene_out))
