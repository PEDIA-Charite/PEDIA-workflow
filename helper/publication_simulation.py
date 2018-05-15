#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import csv 
import getopt
import sys
import argparse
import random
from shutil import copyfile
import glob
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def parse_syn_case_list(file_name):
    # Parse case list
    # 1: omim 2: syndrome, 4: suitable 11: passed ID
    # (0, 'assigned team member')
    # (1, 'Omim ID')
    # (2, 'Syndrome Name')
    # (3, 'F2G: N° of cases')
    # (4, 'suitable for exome analysis')
    # (5, '')
    # (6, 'PEDIA: N° of cases with selected syndrome')
    # (7, 'comment')
    # (8, 'IDs (internal)')
    # (9, 'IDs (external)')
    # (10, 'N° pass QC')
    # (11, 'IDs pass QC')
    # (12, 'internal cases fail QC')
    # 0 cases needed to be newly annotated (min)                                     
    # 1 cases needed
    # 2 assigned to
    # 3 Omim ID
    # 4 Syndrome Name
    # 5 F2G: N° of cases
    # 6 suitable for exome analysis
    # 7 N° of cases pass
    # 8 Cases pass
    # 9 N° of cases fail
    # 10 Cases fail
    # 11 N° of internal cases fail
    # 12 Internal cases fail
    # 13 N° of cases in PEDIA
    # 14 Cases in PEDIA
    # 15 comment

    # init 
    syn_init = {"cases": [], "syndrome_name": "", "omim_id": [], "pub_number": 0, "pedia_number": 0}
    total_num = 0
    suitable_num = 0
    syn_case_dict = {}

    # Parse syndroms and our case list
    with open(file_name, newline='') as csvfile:
        list_reader = csv.reader(csvfile, delimiter='\t')
        for row in list_reader:
            if 'yes' in row[6]:
                syn_name = row[4].replace('\n', ' ')
                cases = row[8].replace('\n', '').split(', ')
                pass_cases = list(filter(lambda x: x!= '0', cases))
                num = int(row[7]) if row[7] != '' else 0
                pub_num = int(row[5]) if row[5] != '' else 0
                suitable_num = suitable_num + num
                syn_case_dict.update({
                    syn_name: {
                        "cases": pass_cases,
                        "syndrome_name": syn_name,
                        "omim_id": row[3],
                        "f2g_number": pub_num,
                        "pedia_number": num,
                        "f2g_id": []
                        }
                    }
                    )

                total_num = total_num + len(pass_cases)
            else:
                syn_name = row[4].replace('\n', ' ')
                no_syn_list.append(syn_name)
        print(suitable_num)
    return syn_case_dict

def parse_f2g_list(file_name, syn_case_dict):
    pub_file = open(file_name)
    pub_data_dict = {}
    pub_data = json.load(pub_file)
    for data in pub_data:
        case_id = data['case_id']
        selected_syn = data['selected_syndromes'][0]
        names = selected_syn['syndrome_name']
        if len(names) == 1:
            name = names[0]
            if name in syn_case_dict:
                syn_case = syn_case_dict[name]
                syn_case['f2g_id'].append(case_id)
                pub_data_dict.update({case_id: data['suggested_syndrome']})
            else:
                if name not in no_syn_list:
                    case_id = data['case_id']
                    print(case_id)
                    print(name)
        else:
            flag = False
            for name in selected_syn['syndrome_name']:
                if name in syn_case_dict:
                    syn_case = syn_case_dict[name]
                    syn_case['f2g_id'].append(case_id)
                    pub_data_dict.update({case_id: data['suggested_syndrome']})
                    flag = True
                    break
                if not flag and name in no_syn_list:
                    flag = True

            if not flag:
                case_id = data['case_id']
                print(case_id)
                print(names)
    return syn_case_dict, pub_data_dict 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate simulated file')
    parser.add_argument('-a', '--aws', help='path to original file')
    parser.add_argument('-p', '--pub', help='path to f2g publication file')
    parser.add_argument('-l', '--list', help='path to PEDIA list file')
    parser.add_argument('-o', '--out', help='path to output')
    parser.add_argument('-c', '--config', help='path to config.yml')

    args = parser.parse_args()
    case_path = os.path.join(args.aws, 'cases')
    pub_file = args.pub
    list_file = args.list
    out_path = args.out
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    log_file = out_path + '/simulation.log'
    logging.basicConfig(filename=log_file, format='%(asctime)s: %(name)s - %(message)s', datefmt='%m-%d %H:%M', level=logging.DEBUG, filemode='w')

    console_handle = logging.StreamHandler()
    console_handle.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(message)s', datefmt='%m-%d %H:%M')
    console_handle.setFormatter(formatter)
    logger.addHandler(console_handle)

    REP_NUM = 10

    logger.debug("Command: %s", str(args))
    logger.info("AWS directory: %s", case_path)
    logger.info("F2G publication file: %s", pub_file)
    logger.info("Deep gestalt list file: %s", list_file)
    logger.info("Output directory: %s", out_path)


    config_file = open(args.config, 'r')
    reader = csv.reader(config_file)
    flag = False
    passed_cases = []
    for row in reader:
        if row[0] == 'TEST_SAMPLES:':
            flag = False
        if row[0] == 'SINGLE_SAMPLES:' or row[0] == 'VCF_SAMPLES:':
            flag = True
            continue
        if flag:
            passed_cases.append(row[0].replace(' ', '').split('-')[1])
    no_syn_list = []
    syn_case_dict = parse_syn_case_list(list_file) 
    syn_case_dict, pub_data_dict = parse_f2g_list(pub_file, syn_case_dict)
    count = 0
    #for syn in syn_case_dict:
    #    print(syn)
    #    print(syn_case_dict[syn])
    #    #count = count + len(syn_case_dict[syn]['f2g_id'])
    #    #if len(syn_case_dict[syn]['cases']) != syn_case_dict[syn]['pedia_number']:
    #    #    print(syn)
    #    #    print(syn_case_dict[syn])
    #    if syn_case_dict[syn]['pedia_number'] > 0:
    #        count = count + len(syn_case_dict[syn]['f2g_id'])
    #print(len(syn_case_dict))
    #            

    origin_vcf_path = os.path.join(args.aws, 'vcfs')
    vcf_list = os.listdir(origin_vcf_path)

    for rep in range(0, REP_NUM):
        logger.info("Start repetition %s", str(rep))
        out_path_rep = os.path.join(out_path, 'REP_' + str(rep))
        if not os.path.exists(out_path_rep):
            os.makedirs(out_path_rep)
        out_case_path = os.path.join(out_path_rep, 'cases')
        if not os.path.exists(out_case_path):
            os.makedirs(out_case_path)
        out_genomic_path = os.path.join(out_path_rep, 'genomics_entries')
        if not os.path.exists(out_genomic_path):
            os.makedirs(out_genomic_path)
        out_vcf_path = os.path.join(out_path_rep, 'vcfs')
        if not os.path.exists(out_vcf_path):
            os.makedirs(out_vcf_path)
        mutation_path = 'data/PEDIA/mutations'
        out_mutation_path = os.path.join(out_path_rep, mutation_path)
        if not os.path.exists(out_mutation_path):
            os.makedirs(out_mutation_path)
        out_simulated_path = os.path.join(out_path_rep, 'simulated_case')
        if not os.path.exists(out_simulated_path):
            os.makedirs(out_simulated_path)

        yml_file = open(os.path.join(out_path_rep, 'config.yml'), "w") 
        yml_file.write("SIMULATE_SAMPLES:\n")
        csvfile = open(os.path.join(out_path_rep, 'simulate_overview.csv'), 'w')
        file_writer = csv.writer(csvfile, delimiter='\t')

        for name in syn_case_dict.keys():
            syn = syn_case_dict[name]
            logger.info("Syndrome name %s", name)
            if syn['pedia_number'] == 0:
                logger.info("%s do not have pedia cases", name)
                continue
            f2g_num = syn['f2g_number']
            case_list = [id for id in syn['cases'] if id in passed_cases]
            failed_case_list = [id for id in syn['cases'] if id not in passed_cases]
            logger.info("Not in passed list: %s", str(failed_case_list))
            pedia_num = len(case_list)

            logger.debug("Number of PEDIA case: %s", str(pedia_num))
            logger.debug("Number of F2G case: %s", str(f2g_num))
            logger.debug("PEDIA case: %s", str(case_list))
            logger.debug("F2G case: %s", str(syn['f2g_id']))
           
            if pedia_num >= f2g_num:
                simulated_case_id = random.sample(case_list, f2g_num)
            else:
                simulated_case_id = [random.choice(case_list) for _ in range(f2g_num)]
            simulated_f2g_id = syn['f2g_id']

            logger.debug('PEDIA cases: %s', str(simulated_case_id))
            logger.debug('F2G cases: %s', str(simulated_f2g_id))

            file_writer.writerow([name, f2g_num, pedia_num, simulated_case_id, simulated_f2g_id, case_list, syn['f2g_id']])

            for case_id, f2g_id in zip(simulated_case_id, simulated_f2g_id):
                file_content = json.load(open(os.path.join(case_path, case_id + '.json')))
                f2g_content = pub_data_dict[f2g_id]

                # Append gestalt to file_content then output
                # Reset gestalt score
                detected_syn = []
                for syndrome in file_content['detected_syndromes']:
                    syndrome['gestalt_score'] = 0
                    detected_syn.append(syndrome)

                count = 0
                for f2g_syn in f2g_content:
                    found = False
                    for syndrome in detected_syn:
                        if syndrome['syndrome_name'] == f2g_syn['syndrome_name']:
                            syndrome['gestalt_score'] = f2g_syn['gestalt_score']
                            found = True
                    if not found:
                        new_syn = f2g_syn
                        new_syn.update({'combined_score': 0, 'feature_score': 0})
                        detected_syn.append(new_syn)

                file_content['detected_syndromes'] = detected_syn
                file_content['f2g_pub_case'] = f2g_id
                file_content['pedia_case_id'] = case_id
                file_content['case_id'] = str(case_id) + '_' + str(f2g_id)

                logger.debug('Simulate: F2G - %s with PEDIA %s', str(f2g_id), case_id)

                with open(os.path.join(out_case_path, str(case_id) + '_' + str(f2g_id) + '.json'), 'w') as outfile:
                    json.dump(file_content, outfile)
                for entries in file_content['genomic_entries']:
                    copyfile(os.path.join(args.aws, 'genomics_entries', str(entries) + '.json'), os.path.join(out_genomic_path, str(entries) + '.json')) 

                for vcf_name in vcf_list:
                    if case_id in vcf_name:
                        copyfile(os.path.join(args.aws, 'vcfs', vcf_name), os.path.join(out_vcf_path, str(case_id) + '_' + str(f2g_id) + '.' + '.'.join(vcf_name.split('.')[1:])))

                cmd = ("python3 preprocess.py -s " + os.path.join(out_case_path, str(case_id) + '_' + str(f2g_id) + '.json') + " -o " + out_simulated_path + ' --skip-vcf')
                logger.debug("Command: %s", cmd)
                os.system(cmd)
                yml_file.write(" - '" + str(case_id) + "_" + str(f2g_id) + "'\n")

                mutation_list = os.listdir(mutation_path)

                #if str(case_id) + '.vcf.gz' in mutation_list:
                #    copyfile(os.path.join(mutation_path, str(case_id) + '.vcf.gz'), os.path.join(out_mutation_path, str(case_id) + '.vcf.gz'))
                #else:
                #    logger.info('%s mutation not found!', str(case_id))

        csvfile.close()
        yml_file.close()
