#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#############################################################################3
# Compare syndrome between f2g and the file we collect
##############################################################################
import json
import os
import csv 
import getopt
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-o', '--output', help='path to output file')
    f2g = '216_gestalt_syn_to_omim.csv'
    pedia = 'test/gestalt_syn_to_omim.csv'
    args = parser.parse_args()
    # syn_dict = {'name':[omim_id]}
    f2g_syn_dict = {}
    pedia_syn_dict = {}
    f2g_file = open(f2g)
    reader = csv.reader(f2g_file, delimiter='\t')
    for row in reader:
        f2g_syn_dict[row[0]] = row[1]
    pedia_file = open(pedia)
    pedia_reader = csv.reader(pedia_file, delimiter='\t')
    for row in pedia_reader:
        pedia_syn_dict[row[0]] = row[1]
    out = {}
    not_out = {}
    for syn in f2g_syn_dict:
        if syn in pedia_syn_dict:
            out[syn] = [pedia_syn_dict[syn], f2g_syn_dict[syn]]
        else:
            not_out[syn] = f2g_syn_dict[syn]
    print(len(out))
            
    syn_file = open(os.path.join(args.output,'same.csv'), 'w')
    syn_writer = csv.writer(syn_file, delimiter='\t')
    for syn in out:
        syn_writer.writerow([syn, out[syn][0], out[syn][1]])
    syn_file.close()
    syn_file = open(os.path.join(args.output,'no_match.csv'), 'w')
    syn_writer = csv.writer(syn_file, delimiter='\t')
    for syn in not_out:
        syn_writer.writerow([syn, not_out[syn]])
    syn_file.close()
