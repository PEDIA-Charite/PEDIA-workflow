# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 12:00:22 2017

@author: Tori
"""
import datetime as dt # Abspeicherung des Files nach Datum ( um Historie verfolgen zu können, eventuell noch Änderung)

import json # JSON öffnen, bearbeiten und speichern 7 open, change and save JSONs
import os, shutil, re # Directory-Informationen bekommen / get information of directory

# change HGVS with mutalyzer
import re

# VCF-Tabelle einrichten / create table for multi-VCF
import csv
import pandas as pd

# GetOpt to read cli inputs
import sys, getopt
import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import hgvs.assemblymapper
import hgvs.exceptions
import requests
import getopt
import sys
import argparse


parser = argparse.ArgumentParser(description='Quality check and generate multi VCF')
parser.add_argument('-m', '--mappedjsons', help='path of mapped json folder')
parser.add_argument('-f', '--vcf', help='path of output multi VCF file')
parser.add_argument('-o', '--output', help='path of currated json files')
parser.add_argument('-s', '--sample', help='path of output sample yml file')
parser.add_argument('-e', '--errorlog', help='path of error log file')
args = parser.parse_args()

jsonoriginal = args.mappedjsons
jsoncurratedfolder = args.output 
mVCF= args.vcf
sample_file = args.sample
err_file = args.errorlog

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37')

file_list = []

multivcf=pd.DataFrame(columns=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', 'NM'])
vcfcounter=0
x=0
var_type_array = []
hgvs_error = open(err_file, 'w')
hgvs_writer = csv.writer(hgvs_error, delimiter='\t')
hgvs_writer.writerow(['id','message'])
vcf_list = []
non_vcf_list = []
test_vcf_list = []
file_counter = 0
seq_list = {}
for filename in os.listdir(jsonoriginal):
    with open(jsonoriginal + '/' + filename) as json_data:
        print(' - ' + filename.split('.json')[0])
        file_counter = file_counter + 1
        d = json.load(json_data)

        gene_list = d['geneList']
        gestalt_found = False
        for gene_entry in gene_list:
            if 'gestalt_score' in gene_entry:
                g_score = gene_entry['gestalt_score']
                if g_score > 0:
                    gestalt_found = True
            if gestalt_found:
                break
        if gestalt_found == False:
            continue
        vcf_found = False
        copy_flag = False
        caseID = d['case_id']


        # This case has a long deletion. It can't be annotated by Jannovar
        if str(caseID) == '158496':
            continue


        if 'original_filename' in d['vcf']:
            vcf_found = True
        hgvslist = []
        if len(d['genomicData']) == 0:
            continue
        for mutation in d['genomicData']:
            if 'Test Information' in mutation:
                if 'Gene Name' in mutation['Test Information']:
                    gene = mutation['Test Information']['Gene Name']
                    if gene == '':
                        if vcf_found == True and caseID not in test_vcf_list:
                            test_vcf_list.append(caseID)
                            copy_flag = True
            if 'HGVS-code' in mutation['Mutations'].keys():
                hgvs=mutation['Mutations']['HGVS-code']
                if hgvs != "":
                    if ", " in hgvs:
                        hgvs_array = hgvs.split(", ")
                        for hgvs in hgvs_array:
                            if hgvs not in hgvslist:
                                hgvslist.append(hgvs)
                    else:
                        hgvslist.append(hgvs)
            elif 'Mutation 1' in mutation['Mutations'].keys():
                for mutationnr, mutationdict in mutation['Mutations'].items():
                    if 'HGVS-code' in mutationdict.keys():
                        hgvs=mutationdict['HGVS-code']
                        hgvslist.append(hgvs)
            else:
                print('keine Mutationen: ', filename)
                if vcf_found == True and caseID not in test_vcf_list:
                    test_vcf_list.append(caseID)
                    copy_flag = True
                continue
            genotype=mutation['Test Information']['Genotype']
            if genotype=='HEMIZYGOUS':
                genotype='1'
            elif genotype=='HOMOZYGOUS':
                genotype='1/1'
            elif genotype=='HETEROZYGOUS' or genotype=='COMPOUND_HETEROZYGOUS':
                genotype='0/1'
            else:
                genotype='0/1'
            if len(hgvslist) == 0:
                if vcf_found == True and caseID not in test_vcf_list:
                    test_vcf_list.append(caseID)
                    copy_flag = True
            for hgvscode in hgvslist:
                try:
                    # Parse HGVS
                    #chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(str(hgvscode), genome, get_transcript=get_transcript)
                    #dup, del, sub, ins
                    variant = hp.parse_hgvs_variant(hgvscode.encode('utf8'))
                    variant_g = vm.c_to_g(variant)
                    variant_g.ac.split('.')[0]
                    chrom = str(int(variant_g.ac.split('.')[0][-2:]))
                    offset = variant_g.posedit.pos.start.base
                    var_type = variant_g.posedit.edit.type
                    found = True
                    if var_type not in var_type_array:
                        var_type_array.append(var_type)
                    if var_type == "identity":
                        print(filename)
                    #if var_type != 'sub':
                       #seq = hdp.get_seq(variant_g.ac)
                    if var_type == 'ins':
                        seq = ""
                        if variant_g.ac in seq_list:
                            seq = seq_list[variant_g.ac]
                        else:
                            seq = hdp.get_seq(variant_g.ac)
                            seq_list.update({variant_g.ac:seq})
                        ref = seq[offset - 1]
                        alt = ref + variant_g.posedit.edit.alt
                    elif var_type == 'dup':
                        seq = ""
                        if variant_g.ac in seq_list:
                            seq = seq_list[variant_g.ac]
                        else:
                            seq = hdp.get_seq(variant_g.ac)
                            seq_list.update({variant_g.ac:seq})
                        end = variant_g.posedit.pos.end.base
                        offset = end
                        ref = seq[end - 1]
                        alt = ref + variant_g.posedit.edit.ref
                    elif var_type == 'del':
                        seq = ""
                        if variant_g.ac in seq_list:
                            seq = seq_list[variant_g.ac]
                        else:
                            seq = hdp.get_seq(variant_g.ac)
                            seq_list.update({variant_g.ac:seq})

                        offset = offset - 1
                        ref = variant_g.posedit.edit.ref
                        ref = seq[offset - 1] + ref
                        alt = seq[offset - 1]
                    elif var_type == 'sub' or var_type == 'delins':
                        ref = variant_g.posedit.edit.ref
                        alt = variant_g.posedit.edit.alt
                    else:
                        print(hgvscode)
                        print(var_type)
                        found = False
                        
                    if found == True:
                        if vcf_found == False and caseID not in non_vcf_list:
                            non_vcf_list.append(caseID)
                            copy_flag = True
                        if vcf_found == True and caseID not in vcf_list:
                            vcf_list.append(caseID)
                            copy_flag = True
                        if hgvscode in multivcf['NM'].tolist():
                            index=multivcf['NM'].tolist().index(hgvscode)
                            multivcf.set_value(index, caseID, genotype)
                        else:
                            #chromo=chrom.split('chr')[1]
                            chromo=chrom
                            print(type(chrom))
                            if chrom == '23':
                                chrom = 'X'
                            if chrom == '24':
                                chrom = 'Y'
                            multivcf.set_value(x, '#CHROM', str(chrom))
                            try:
                                multivcf.set_value(x,'sort',int(chromo))
                            except ValueError as e:
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

                except Exception as e:
                    pass
                    hgvs_writer.writerow([caseID, str(e)])
        if copy_flag == True:
            cmd = 'cp -f ' + jsonoriginal + '/' + filename + ' ' + jsoncurratedfolder
            os.system(cmd)
print(var_type_array)
hgvs_error.close()

##data_vcf sortieren
print('Sort DataFrame ...')
multivcf=multivcf.sort_values(by=['sort', 'POS'])
multivcf=multivcf.drop('sort',axis=1)
multivcf=multivcf.drop('NM',axis=1)
multivcf=multivcf.reset_index(drop=True)

#leere Felder füllen
multivcf=multivcf.fillna(value='0/0')

multivcf.to_csv(mVCF+".tmp", sep='\t', index=False, header=True, quoting=csv.QUOTE_NONE)

with open(mVCF, 'w') as outfile:
    outfile.write('##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    outfile.write('##contig=<ID=1,assembly=b37,length=249250621>\n')
    outfile.write('##contig=<ID=2,assembly=b37,length=243199373>\n')
    outfile.write('##contig=<ID=3,assembly=b37,length=198022430>\n')
    outfile.write('##contig=<ID=4,assembly=b37,length=191154276>\n')
    outfile.write('##contig=<ID=5,assembly=b37,length=180915260>\n')
    outfile.write('##contig=<ID=6,assembly=b37,length=171115067>\n')
    outfile.write('##contig=<ID=7,assembly=b37,length=159138663>\n')
    outfile.write('##contig=<ID=8,assembly=b37,length=146364022>\n')
    outfile.write('##contig=<ID=9,assembly=b37,length=141213431>\n')
    outfile.write('##contig=<ID=10,assembly=b37,length=135534747>\n')
    outfile.write('##contig=<ID=11,assembly=b37,length=135006516>\n')
    outfile.write('##contig=<ID=12,assembly=b37,length=133851895>\n')
    outfile.write('##contig=<ID=13,assembly=b37,length=115169878>\n')
    outfile.write('##contig=<ID=14,assembly=b37,length=107349540>\n')
    outfile.write('##contig=<ID=15,assembly=b37,length=102531392>\n')
    outfile.write('##contig=<ID=16,assembly=b37,length=90354753>\n')
    outfile.write('##contig=<ID=17,assembly=b37,length=81195210>\n')
    outfile.write('##contig=<ID=18,assembly=b37,length=78077248>\n')
    outfile.write('##contig=<ID=19,assembly=b37,length=59128983>\n')
    outfile.write('##contig=<ID=20,assembly=b37,length=63025520>\n')
    outfile.write('##contig=<ID=21,assembly=b37,length=48129895>\n')
    outfile.write('##contig=<ID=22,assembly=b37,length=51304566>\n')
    outfile.write('##contig=<ID=X,assembly=b37,length=155270560>\n')
    outfile.write('##contig=<ID=Y,assembly=b37,length=59373566>\n')
    with open(mVCF+".tmp",'r') as infile:
        for line in infile:
            outfile.write(line)

os.remove(mVCF+".tmp")

with open(sample_file, 'w') as sample_out:
    sample_out.write('SINGLE_SAMPLES:\n')
    for case in non_vcf_list:
        sample_out.write(' - ' + str(case) + '\n')
    sample_out.write('VCF_SAMPLES:\n')
    for case in vcf_list:
        sample_out.write(' - ' + str(case) + '\n')
    sample_out.write('TEST_SAMPLES:\n')
    for case in test_vcf_list:
        sample_out.write(' - ' + str(case) + '\n')
