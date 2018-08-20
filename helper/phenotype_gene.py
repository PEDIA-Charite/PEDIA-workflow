#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import csv 
import getopt
import sys
import argparse

def get_case_info(file_writer, file_content, reason, status="Pass"):
    out = []
    tmp_out = []
    info = []
    ge_out = []
    hgvs_out = []
    final_syn = []
    all_gene = []
    all_omim = []
    status = ''
    syndrome_name = "" 
    syndrome = {} 
    has_mask = ''
    if 'genomicData' in file_content:
        genomic = file_content['genomicData']
        ge_out = [ge_entry['Test Information']['Gene Name'] for ge_entry in genomic]


    # Check selected syndrome
    if 'selected_syndromes' in file_content:
        syn_out = []
        omim_out = []
        diag_out = []
        has_mask_out = []
        for syn in file_content['selected_syndromes']:
            if syn['diagnosis'] == 'DIFFERENTIAL_DIAGNOSIS':
                continue
            syndrome.update({syn['syndrome_name']: {
                'omim': syn['omim_id'],
                'diagnosis': syn['diagnosis'],
                'has_mask': syn['has_mask']
                }})
            syn_out.append(syn['syndrome_name'])
            omim_out.append(syn['omim_id'])
            diag_out.append(syn['diagnosis'])
            has_mask_out.append(syn['has_mask'])

        if len(syn_out) == 1:
            # only one syndrome selected
            syndrome_name = syn_out[0]
            has_mask = has_mask_out[0]
            final_syn.append(omim_out[0])

            # Store omim id into all_omim
            if type(omim_out[0]) == list:
                for omim in omim_out[0]:
                    all_omim.append(omim)
            else:
                all_omim.append(omim_out[0])
        else:
            list_syn = []
            single_syn = []
            for idx, syn in enumerate(omim_out):
                if type(syn) == list:
                    list_syn = syn
                    syndrome_name = syn_out[idx]
                else:
                    single_syn.append(syn)

            final_syn.append(list_syn)
            all_omim = [omim for omim in list_syn]
            for syn in single_syn:
                if syn not in list_syn:
                    final_syn.append(syn)
                    all_omim.append(syn)

    # Fixed deprecated omim id
    tmp_all_omim = []
    for omim in all_omim: 
        if str(omim) in deprecated_omim:
            for new_id in deprecated_omim[str(omim)]:
                tmp_all_omim.append(new_id)
        else:
            tmp_all_omim.append(omim)
    all_omim = tmp_all_omim
    # Map phenotype to gene
    for omim in all_omim:
        if str(omim) in omim_dict:
            genes = omim_dict[str(omim)]
            for gene in genes:
                if gene not in all_gene:
                    all_gene.append(gene)

    # Syndrome has no omim id
    if len(all_omim) == 0:
        all_gene = [gene for gene in ge_out]
        status = 'no omim'
    # no gene can be found via mapping 
    if len(all_gene) == 0:
        all_gene = [gene for gene in ge_out]
        status = 'no gene mapped from omim'
    if ge_out[0] not in all_gene:
        status = 'disease gene not in omim gene list'

    if syndrome_name not in tmp_syn:
        syn_cases = []
        cases_count = 0
    else:
        syn_cases = tmp_syn[syndrome_name]['cases']
        cases_count = tmp_syn[syndrome_name]['count']
    syn_cases.append(fileName.split('.')[0])
    cases_count = cases_count + 1
    tmp_syn.update({syndrome_name:{
        'gene':all_gene,
        'omim':all_omim,
        'has_mask':syndrome[syndrome_name]['has_mask'],
        'cases': syn_cases,
        'count': cases_count
        }})

    # Output results to csv
    tmp_out.append(fileName.split('.')[0])
    tmp_out.append(status)
    tmp_out.append(has_mask_out)
    tmp_out.append(1 if 1 in has_mask_out else 0)
    tmp_out.append(all_gene)
    tmp_out.append(str(ge_out)[1:-1])
    tmp_out.append(final_syn)
    tmp_out.append(syn_out)
    tmp_out.append(diag_out)
    tmp_out.append(omim_out)
    file_writer.writerow(tmp_out)
    if status != '':
        not_mapped_writer.writerow(tmp_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate summary file')
    parser.add_argument('-c', '--case', help='path to convert file')
    parser.add_argument('-l', '--log', help='path to log file')
    parser.add_argument('-o', '--output', help='path to output file')

    args = parser.parse_args()
    case_path = args.case
    log_file = open(args.log)
    log_data = json.load(log_file)

    deprecated_file = open('data/omim_deprecated_replacement.json', 'r')
    deprecated_omim = json.load(deprecated_file)
    

    # Parse omim phenotype to gene mapping file
    omim_dict = {} 
    omim_file = open('data/morbidmap.txt', 'r')
    omim_reader = csv.reader(omim_file, delimiter='\t')
    count = 0
    for row in omim_reader:
        count = count + 1
        if len(row) > 2:
            omim = row[0].split(', ')[-1].split(' ')[0]
            gene = row[1].split(', ')[0]
            if omim not in omim_dict:
                omim_dict.update({omim:[gene]})
            else:
                if gene not in omim_dict[omim]:
                    omim_dict[omim].append(gene)

    omim_file.close()
    output_filename = os.path.join(args.output, 'case_gene_phenotype_table.csv')

    not_mapped_file = open(os.path.join(args.output,'gene_not_mapped.csv'), 'w')
    not_mapped_writer = csv.writer(not_mapped_file, delimiter='\t')
    not_mapped_writer.writerow(['case', 'status', 'has_mask_all_syndrome', 'has_mask', 'all genes', 'gene', 'omim', 'syndrome name', 'diagnosis', 'all omim', 'count'])

    syn_mapped_file = open(os.path.join(args.output,'phenotype_to_gene.csv'), 'w')
    syn_mapped_writer = csv.writer(syn_mapped_file, delimiter='\t')
    tmp_syn = {}

    with open(output_filename, 'w') as csvfile:
        file_writer = csv.writer(csvfile, delimiter='\t')
        file_writer.writerow(['case', 'status','has_mask_all_syndrome', 'has_mask', 'all genes', 'gene', 'omim', 'syndrome name', 'diagnosis', 'all omim'])

        for fileName in log_data["passed"].keys():
            test = fileName
            file_content = json.load(open(os.path.join(case_path, fileName + '.json')))
            get_case_info(file_writer, file_content, [])

    final_syn = {}
    for syn in tmp_syn:
        if len(tmp_syn[syn]['omim']) > 1 or len(tmp_syn[syn]['omim']) == 0:
            final_syn.update({syn:tmp_syn[syn]})

    for syn in tmp_syn:
        if len(tmp_syn[syn]['omim']) == 1:
            found = False
            for syn_compare in tmp_syn:
                if len(tmp_syn[syn_compare]['omim']) > 1:
                    if tmp_syn[syn]['omim'][0] in tmp_syn[syn_compare]['omim']:
                        for case in tmp_syn[syn]['cases']: 
                            final_syn[syn_compare]['cases'].append(case)
                        final_syn[syn_compare]['count'] = final_syn[syn_compare]['count'] + tmp_syn[syn]['count']
                        found = True
                        break
            if not found:
                final_syn.update({syn:tmp_syn[syn]})

    all_genes = []
    syn_mapped_writer.writerow(['name', 'has_mask', 'gene', 'omim', 'cases', 'count'])
    for syn in final_syn:
        syn_mapped_writer.writerow([syn, tmp_syn[syn]['has_mask'], tmp_syn[syn]['gene'], tmp_syn[syn]['omim'], tmp_syn[syn]['cases'], tmp_syn[syn]['count']])
        for gene in tmp_syn[syn]['gene']:
            if gene not in all_genes:
                all_genes.append(gene)

    gene_mapped_file = open(os.path.join(args.output,'gene_mapped_table.csv'), 'w')
    gene_mapped_writer = csv.writer(gene_mapped_file, delimiter='\t')
    for gene in all_genes:
        gene_mapped_writer.writerow([gene])

    gene_mapped_file.close()
    not_mapped_file.close()
    syn_mapped_file.close()
