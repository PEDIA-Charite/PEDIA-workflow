#!/bin/env python3
import json
import sys

with open(sys.argv[1], "r") as jsfile:
    data = json.load(jsfile)

for e in data["genomic_entries"]:
    print(e)
