#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: hertzberg,leitheim
"""

import json
import os
import csv  # necessary for creating genedict
import requests
import getopt
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mapping disorder to gene')
    parser.add_argument('-j', '--jsonsoriginal', help='path of original json folder')
    parser.add_argument('-u', '--usi', help='USI number')
    parser.add_argument('-o', '--output', help='path of mapped json folder')
    args = parser.parse_args()

    path = args.jsonsoriginal
    usi = args.usi
    file_content = json.load(open(path))
    file_content['vcf'] = [usi + '.vcf.gz']
    out_file = os.path.join(args.output, path.split('/')[-1])
    with open(out_file, 'w') as outfile:
        json.dump(file_content, outfile)
