import os
import json
import csv
from pprint import pprint
from collections import defaultdict

with open("f2g_library_dump.json", "r") as f2g_lib:
    syndromes = json.load(f2g_lib)

ps_dict = {}
mim_to_ps = {}
with open("data/phenotypicSeries.txt", "r") as ps_data:
    for raw_entry in csv.DictReader(
            (r for r in ps_data if not r.startswith("#")),
            delimiter="\t",
            fieldnames=["ps_number", "mim_number", "phenotype"]
    ):
        if raw_entry["ps_number"] not in ps_dict:
            ps_dict[raw_entry["ps_number"]] = {
                "name": raw_entry["mim_number"],
                "syndromes": []
            }
        else:
            ps_dict[raw_entry["ps_number"]]["syndromes"].append(
                {
                    "mim_number": raw_entry["mim_number"],
                    "name": raw_entry["phenotype"]
                }
            )
            if raw_entry["mim_number"] in mim_to_ps:
                mim_to_ps[raw_entry["mim_number"]].append(
                    raw_entry["ps_number"]
                )
            else:
                mim_to_ps[raw_entry["mim_number"]] = [raw_entry["ps_number"]]

f2g_only_series = [
    "PS104300",
    "PS143465",
    "PS612863",
    "PS121700",
    "PS125853",
    "PS601551",
    "PS164210",
    "PS235000",
    "PS145500",
    "PS608516",
    "PS166710",
]

f2g_names = {}

name_to_sid = {}

errors = []

name_to_entry = {}

json_files = [v for v in os.listdir("process/aws_dir/cases")
              if v.endswith(".json")]

for syndrome in syndromes:
    sid = syndrome["syndrome_id"]
    oid = syndrome["omim_id"]
    psid = syndrome["omim_ps_id"]
    if psid and not psid.startswith("PS"):
        psid = "PS" + psid
    sname = syndrome["syndrome_name"]
    salt = [s["synonym"] for s in syndrome["synonyms"]]
    salt += [sname]

    nbk_id = syndrome["nbk_id"]
    orphanet = syndrome["orphanet"]

    if sname not in name_to_entry:
        name_to_entry[sname.lower()] = syndrome
    # else:
    #     print(sname)
    #     print(name_to_entry[sname])
    #     print(syndrome)

    if oid == "0" or oid is None:
        # should be names of phenotypic series
        # if psid is None:
        #     print(oid, sname, nbk_id, orphanet)
        continue

trans = {
    "614132": "213980"
}

def match_syndrome(name: str, oid) -> str:
    try:
        sname = name.lower()
    except AttributeError as e:
        print(name, e)

    oid = [oid] if not isinstance(oid, list) else oid

    oid = [str(o) for o in oid]
    oid = [trans[o] if o in trans else o for o in oid]

    if len(oid) < 2:
        return None

    if sname in name_to_entry:
        return name_to_entry[sname]

    splits = sname.split(";")
    for s in splits:
        s = s.strip()
        for k in name_to_entry:
            if s in k:
                soid = name_to_entry[k]["omim_id"]
                if soid in oid:
                    return name_to_entry[k]
                else:
                    print("{} not in {}".format(soid, oid))

    print(sname, "no match")

matched = 0
for js_file in json_files:
    js_path = "process/aws_dir/cases/" + js_file
    with open(js_path, "r") as jfp:
        js_data = json.load(jfp)

    print(js_data["case_id"])

    selected = js_data["selected_syndromes"]
    detected = js_data["detected_syndromes"]
    for d in detected:
        sname = d["syndrome_name"]
        oid = d["omim_id"]
        r = match_syndrome(sname, oid)
        if not None:
            matched += 1

print("Matched", matched)
