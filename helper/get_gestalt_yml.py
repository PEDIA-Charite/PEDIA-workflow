#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import csv 
import getopt
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate config.yml file with gestalt scores')
    parser.add_argument('-i', '--input', help='path to phenotype_to_gene.csv')
    parser.add_argument('-o', '--output', help='path to output file')

    args = parser.parse_args()
    input_path = args.input

    # Parse omim phenotype to gene mapping file
    input_file = open(input_path, 'r')
    input_reader = csv.reader(input_file, delimiter='\t')
    cases = []
    count = 0
    for row in input_reader:
        count = count + 1
        if row[1] == '1':
            tmp = row[4][1:-1].split(', ')
            print(tmp)
            for case in tmp:
                if case not in cases:
                    cases.append(case)
    output_filename = args.output
    with open(output_filename, 'w') as csvfile:
        csvfile.write('GESTALT_SAMPLES:\n')
        #file_writer = csv.writer(csvfile, delimiter='\t')
        for case in cases:
            csvfile.write('- ' + case[1:-1] + '\n')


