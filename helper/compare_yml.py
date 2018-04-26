import csv
import os 
import json
import sys
import argparse

def diff(first, second):
    second = set(second)
    return [item for item in first if item not in second]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='compare config')
    parser.add_argument('-p', '--previous', help='path of original config.yml')
    parser.add_argument('-n', '--new', help='path of new config.yml')
    parser.add_argument('-o', '--out', help='path of output')
    args = parser.parse_args()
    old = args.previous
    new = args.new
    out = args.out
    old_file = open(old, 'r')
    new_file = open(new, 'r')
    reader = csv.reader(old_file)
    old_single = []
    new_single = []
    flag = False
    for row in reader:
        if row[0] == 'TEST_SAMPLES:':
            flag = False
        if flag:
            old_single.append(row[0][3:])
        if row[0] == 'SINGLE_SAMPLES:':
            flag = True
    flag = False
    reader2 = csv.reader(new_file)
    for row in reader2:
        if row[0] == "TEST_SAMPLES:":
            flag = False
        if flag:
            new_single.append(row[0][3:])
        if row[0] == "SINGLE_SAMPLES:":
            flag = True
    diff1 = diff(old_single, new_single)
    diff2 = diff(new_single, old_single)
    old_file.close()
    new_file.close()
    outfile = open(out, 'w')
    outfile.write('In old but not in new\n')
    for value in diff1:
        outfile.write(value + '\n')
    outfile.write('In new but not in old\n')
    for value in diff2:
        outfile.write(value + '\n')
    outfile.close()

