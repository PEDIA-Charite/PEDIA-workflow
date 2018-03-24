import json

entries = json.load(open('case_gene_name.json', 'r'))

case_dict = {}
for entry in entries:
    for k, v in entry.items():
        if k in case_dict:
            if v:
                case_dict[k] = v
        else:
            case_dict[k] = v

filtered = [{k: v} for k,v in case_dict.items()]

json.dump(filtered, open('new_names.json', 'w'))
