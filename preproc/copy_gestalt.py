import os
import csv 
import getopt
import sys
import argparse
import shutil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-i', '--case', help='path to convert file')
    parser.add_argument('-o', '--output', help='path to output file')

    args = parser.parse_args()
    case_path = args.case
    name = 'config_gestalt.csv'
    g_file = open(name, 'r')
    reader = csv.reader(g_file)
    file_list = []
    for row in reader:
        file_list.append(row[0] + '.json')

    for case in file_list:
        shutil.copy2(os.path.join(case_path, case), os.path.join(args.output, case))

