# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 01:06:29 2016

@author: Tori
"""

### Multi-VCF aus JSONs mit korrektem HGVS-Code

import json
import os
import csv
import pandas as pd
import pyhgvs
import pyhgvs.utils as hgvs_utils
from pygr.seqdb import SequenceFileDB

## HGVS-Library-Zeugs
# Read genome sequence using pygr.
genome = SequenceFileDB('genome.fa')

# Read RefSeq transcripts into a python dict.
with open('genes.refGene') as infile:
    transcripts = hgvs_utils.read_transcripts(infile)

# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)

## MultiVCF

multivcf=pd.DataFrame(columns=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', 'NM'])

x=0

for file in os.listdir('C:/Users/Tori/Documents/Python Scripts/neue jsons/'):
    if 'json' in file:
        with open('C:/Users/Tori/Documents/Python Scripts/neue jsons/'+file) as json_data:
            d=json.load(json_data)
            caseID=d['case_id']
            hgvslist=[]
            if len(d['genomicData'])!=0:
                for mutation in d['genomicData']:
                    if 'HGVS-code' in mutation['Mutations'].keys():
                        hgvs=mutation['Mutations']['HGVS-code']
                        hgvslist.append(hgvs)
                    elif 'Mutation 1' in mutation['Mutations'].keys():
                        for mutationnr, mutationdict in mutation['Mutations'].items():
                            if 'HGVS-code' in mutationdict.keys():
                                hgvs=mutationdict['HGVS-code']
                                hgvslist.append(hgvs)
                    else:
                        continue
                    genotype=mutation['Test Information']['Genotype']
                    if genotype=='Hemizygous':
                        genotype='1'
                    elif genotype=='Homozygous':
                        genotype='1/1'
                    elif genotype=='Heterozygous' or genotype=='Compound Heterozygous':
                        genotype='0/1'
                    else:
                        genotype='./1'
                    for hgvscode in hgvslist:
                        try:
                            chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
                                str(hgvscode), genome, get_transcript=get_transcript)
                            if hgvscode in multivcf['NM'].tolist():
                                index=multivcf['NM'].tolist().index(hgvscode)
                                multivcf.set_value(index, caseID, genotype)
                            else:
                                chromo=chrom.split('chr')[1]
                                multivcf.set_value(x, '#CHROM', str(chromo))
                                try:
                                    multivcf.set_value(x,'sort',int(chromo))
                                except ValueError,e:
                                    multivcf.set_value(x,'sort',30)
                                multivcf.set_value(x, 'NM', hgvscode)
                                multivcf.set_value(x, 'POS', offset)
                                multivcf.set_value(x, 'ID', '.')
                                multivcf.set_value(x, 'REF', str(ref))
                                multivcf.set_value(x, 'ALT', str(alt))
                                multivcf.set_value(x, 'QUAL', '.')
                                multivcf.set_value(x, 'FILTER', '.')
                                multivcf.set_value(x, 'INFO', 'HGVS="'+hgvscode+'"')
                                multivcf.set_value(x, 'FORMAT', 'GT')
                                multivcf.set_value(x, caseID, genotype)
                                x=x+1
                        except ValueError, e: #'falsche' HGVS-Codes überspringen und anzeigen
                            print 'Error:',file, hgvs, e
                            continue
                        
                        
##data_vcf sortieren
print 'Sort DataFrame ...'
multivcf=multivcf.sort_values(by=['sort', 'POS'])
multivcf=multivcf.drop('sort',axis=1)
multivcf=multivcf.drop('NM',axis=1)
multivcf=multivcf.reset_index(drop=True)

#leere Felder füllen
multivcf=multivcf.fillna(value='0/0')

multivcf.to_csv('mutationsJSONs.vcf', sep='\t', index=False, header=True, quoting=csv.QUOTE_NONE)

with open('JsonsVCF.vcf', 'w') as outfile:
    outfile.write('##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    with open('mutationsJSONs.vcf','r') as infile:
        for line in infile:
            outfile.write(line)            