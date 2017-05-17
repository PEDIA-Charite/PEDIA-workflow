#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, getopt
import json
import glob
import csv
import numpy as np

argv = sys.argv[1:]
	
inputfile = ''
outputfile = ''
try:
	opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
	print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
		sys.exit()
	elif opt in ("-i", "--ifolder"):
		inputfile = arg
	elif opt in ("-o", "--ofile"):
		outputfile = arg
print('Input folder is ',inputfile)
print('Output file is ',outputfile)

table = []

for filename in glob.glob(inputfile+'/*.json'):
	r = open(filename,'r',encoding='ISO-8859-1')
	data = json.loads(r.read())
	case = data['case_id']
	gene = data["genomicData"][0]["Test Information"]["Gene Name"]
	for entry in data['geneList']:
		if entry["gene_symbol"] == gene:
			label = 1
		else:
			label = 0
		table.append({"case": case, "gene_id": entry['gene_id'], "feature_score": entry.get("feature_score", np.nan), "cadd_phred_score": entry.get("cadd_phred_score", np.nan), "combined_score":  entry.get("combined_score", np.nan), "cadd_raw_score":  entry.get("cadd_raw_score", np.nan), "gestalt_score":  entry.get("gestalt_score", np.nan), "boqa_score":  entry.get("boqa_score", np.nan), "pheno_score": entry.get("pheno_score", np.nan), "label": label})
			



with open(outputfile, 'w') as csvfile:
	fieldnames = ["case", "gene_id", "feature_score", "cadd_phred_score", "combined_score", "cadd_raw_score", "gestalt_score", "boqa_score", "pheno_score", "label"]
	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
	writer.writeheader()
	for row in table:
		writer.writerow(row)
