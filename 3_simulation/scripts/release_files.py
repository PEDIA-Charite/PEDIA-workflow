#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import sys
import argparse
import datetime

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mapping disorder to gene')
    parser.add_argument('-i', '--input', help='path of original json folder')
    parser.add_argument('-o', '--output', help='path of output json folder')
    args = parser.parse_args()

    exe_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    version = 'v1.1'

    # Copy CV
    path = os.path.join(args.input, 'jsons/1KG/CV_gestalt')
    out_path = os.path.join(args.output, 'pedia_jsons_v1.1/jsons/1KG/CV_gestalt')
    os.makedirs(out_path, exist_ok=True)
    all_files = os.listdir(path)
    for json_file in all_files:
        file_content = json.load(open(os.path.join(path, json_file)))
        del file_content['submitter']
        del file_content['detected_syndromes']
        out_file = os.path.join(out_path, json_file)
        with open(out_file, 'w') as outfile:
            json.dump(file_content, outfile)

    # Copy real and test
    path = os.path.join(args.input, 'jsons/real/gestalt/train/1KG')
    out_path = os.path.join(args.output, 'pedia_jsons_v1.1/jsons/real/gestalt/train/1KG')
    os.makedirs(out_path, exist_ok=True)
    all_files = os.listdir(path)
    for json_file in all_files:
        file_content = json.load(open(os.path.join(path, json_file)))
        del file_content['submitter']
        del file_content['detected_syndromes']
        out_file = os.path.join(out_path, json_file)
        with open(out_file, 'w') as outfile:
            json.dump(file_content, outfile)

    path = os.path.join(args.input, 'jsons/real/gestalt/test')
    out_path = os.path.join(args.output, 'pedia_jsons_v1.1/jsons/real/gestalt/test')
    os.makedirs(out_path, exist_ok=True)
    all_files = os.listdir(path)
    for json_file in all_files:
        file_content = json.load(open(os.path.join(path, json_file)))
        del file_content['submitter']
        del file_content['detected_syndromes']
        out_file = os.path.join(out_path, json_file)
        with open(out_file, 'w') as outfile:
            json.dump(file_content, outfile)

    log_file = os.path.join(args.output, 'pedia_jsons_v1.1/log')
    with open(log_file, "w") as text_file:
        text_file.write("Time: %s\n" % exe_time)
        text_file.write("Version: %s\n" % version)
        text_file.write("Input foler: %s\n" % args.input)
        text_file.write("Output foler: %s\n" % args.output)

    # Copy publication test
    path = os.path.join(args.input, 'publication_simulation')
    out_path = os.path.join(args.output, 'publication_simulation_v1.1')
    os.makedirs(out_path, exist_ok=True)
    types = ['train', 'CV']
    data = ['1KG']
    out_types = ['train', 'test']
    for i in range(10):
        for j in data:
            for idx, k in enumerate(types):
                sim_path = os.path.join(path, 'REP_' + str(i), 'jsons', j, k)
                all_files = os.listdir(sim_path)
                for json_file in all_files:
                    file_content = json.load(open(os.path.join(sim_path, json_file)))
                    del file_content['detected_syndromes']
                    del file_content['submitter']
                    file_content["algo_deploy_version"] = "280618"
                    out_file_path = os.path.join(out_path, 'REP_' + str(i), 'jsons', j, out_types[idx])
                    os.makedirs(out_file_path, exist_ok=True)
                    out_file = os.path.join(out_file_path, json_file)
                    with open(out_file, 'w') as outfile:
                        json.dump(file_content, outfile)

    log_file = os.path.join(args.output, 'publication_simulation_v1.1/log')
    with open(log_file, "w") as text_file:
        text_file.write("Time: %s\n" % exe_time)
        text_file.write("Version: %s\n" % version)
        text_file.write("Input foler: %s\n" % args.input)
        text_file.write("Output foler: %s\n" % args.output)
