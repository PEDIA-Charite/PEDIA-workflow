#!/bin/env python3
import os
import json
import sys
from pprint import pprint

jspath = sys.argv[1]

case_path, filename = os.path.split(jspath)

base_path, _ = os.path.split(case_path)


with open(jspath, "r") as jsfile:
    data = json.load(jsfile)

for e in data["genomic_entries"]:
    print(e)

    ge_path = os.path.join(base_path, "genomics_entries", "{}.json".format(e))
    if os.path.exists(ge_path):
        with open(ge_path) as gefile:
            ge_data = json.load(gefile)
        pprint(ge_data)
    else:
        print("Could not find entry in path: {}".format(ge_path))
