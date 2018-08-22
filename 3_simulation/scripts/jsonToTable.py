#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, getopt
import json
import glob
import csv
import numpy as np


# Options

argv = sys.argv[1:]

inputfile = ''
outputfile = ''
try:
	opts, args = getopt.getopt(argv,"hi:o:",["help","ifile=","ofile="])
except getopt.GetoptError:
	print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
	sys.exit(2)
for opt, arg in opts:
	if opt in ("-h", "--help"):
		print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
		sys.exit(1)
	elif opt in ("-i", "--ifolder"):
		inputfile = arg
	elif opt in ("-o", "--ofile"):
		outputfile = arg
print('Input folder is ',inputfile)
print('Output file is ',outputfile)


# Load data and use max for duplicates
## key is sampleid_entrezid
def getKey(case,gene_id):
	return str(case) + "_" + str(gene_id)

hashedData = {}

for filename in glob.glob(inputfile+'/*.json'):
	found = False
	geneIDs = []
	r = open(filename,'r',encoding='ISO-8859-1')
	data = json.loads(r.read())
	case = data['case_id']
	gene = data["genomicData"][0]["Test Information"]["Gene Name"]
	for entry in data['geneList']:
		if entry["gene_symbol"] == gene:
			label = 1
			found = True
		else:
			label = 0
		geneID = entry['gene_id']
		geneIDs.append(geneID)
		featureScore = entry.get("feature_score", np.nan)
		caddPhredScore = entry.get("cadd_phred_score", np.nan)
		combinedScore = entry.get("combined_score", np.nan)
		caddRawSscore = entry.get("cadd_raw_score", np.nan)
		gestaltScore = entry.get("gestalt_score", np.nan)
		boqaScore = entry.get("boqa_score", np.nan)
		phenoScore = entry.get("pheno_score", np.nan)


		if getKey(case, geneID) in hashedData:
			value = hashedData[getKey(case, geneID)]
			featureScore = max(featureScore, value["feature_score"])
			caddPhredScore = max(caddPhredScore, value["cadd_phred_score"])
			combinedScore = max(combinedScore, value["combined_score"])
			caddRawSscore = max(caddRawSscore, value["cadd_raw_score"])
			gestaltScore = max(gestaltScore, value["gestalt_score"])
			boqaScore = max(boqaScore, value["boqa_score"])
			phenoScore = max(phenoScore, value["pheno_score"])

		hashedData[getKey(case,geneID)] = {"case": case, "gene_id": geneID, "feature_score": featureScore, "cadd_phred_score": caddPhredScore, "combined_score":  combinedScore, "cadd_raw_score":  caddRawSscore, "gestalt_score":  gestaltScore, "boqa_score":  boqaScore, "pheno_score": phenoScore, "label": label}

	if not found:
		print("problem with",case,gene)
		for geneID in geneIDs:
			if getKey(case,geneID) in hashedData:
				del hashedData[getKey(case,geneID)]

with open(outputfile, 'w') as csvfile:
	fieldnames = ["case", "gene_id", "feature_score", "cadd_phred_score", "combined_score", "cadd_raw_score", "gestalt_score", "boqa_score", "pheno_score", "label"]
	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
	writer.writeheader()
	for key, value in hashedData.items():
		writer.writerow(value)
