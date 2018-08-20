#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv 
import getopt
import sys
import argparse
import shutil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-c', '--case', help='path to convert file')
    parser.add_argument('-o', '--output', help='path to output file')
    parser.add_argument('-l', '--gene-list', help='list of genes you want to find in case')
    
    args = parser.parse_args()
    case_path = args.case
    list_file = open(args.gene_list, 'r')
    reader = csv.reader(list_file, delimiter='\t')
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    for row in reader:
        shutil.copy2(os.path.join(case_path, row[0] + '.json'), os.path.join(args.output, row[0] + '.json'))
    

