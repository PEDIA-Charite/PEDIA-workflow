import json
import os
from pprint import pprint

'''
Discover some structure in json files. Statically programmed. Directly change the script to discover some other specific aspect of the json files. Might be useful in the implementation of json handlers.
'''

dirpath = '/data/users/Max/PEDIA-workflow/process/aws_dir/genomics_entries'

js_files = [ json.load(open(os.path.join(dirpath,o),'r')) for o in os.listdir(dirpath) if o.endswith('.json') ]

keys = {}

count = 0
for j in js_files:
    j = j['variants']
    mutations = []
    if 'mutation' in j:
        mutations.append(j['mutation'])
    elif 'mutation1' in j:
        mutations += [j['mutation1'],j['mutation2']]
    else:
        continue
    count += 1
    for m in mutations:
        if m:
            if 'mutation_type' not in m.keys():
                print(m.keys())
            for k,v in m.items():
                if k in keys:
                    keys[k] += 1
                else:
                    keys[k] = 1

pprint(keys)

print('Count all with mutations', count)

# pprint(abnormal)
