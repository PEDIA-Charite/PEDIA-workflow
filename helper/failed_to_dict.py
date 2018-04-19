import json
import csv

from pprint import pprint

with open("hgvs_errors.json") as hgvs_file:
    hgvs_errors = json.load(hgvs_file)

hgvs_errors["version"] += 1

with open("failedhgvscodes_corrected.txt", encoding="latin1") as new_file:
    dr = csv.DictReader(new_file, delimiter="\t")
    hgvs_new = [{"case_id": c["Case ID"], "corrected": c["corrected HGVS"]}
                for c in dr]

sure_new = {
    str(v["case_id"]): v["corrected"]
    for v in hgvs_new
    if v["corrected"] != '' and "probably" not in v["corrected"]
}

pprint(sure_new)

for s, v in sure_new.items():
    update = {
        "cleaned": [v],
        "note": "failedjannovar_corrected"
    }
    if s in hgvs_errors["data"]:
        hgvs_errors["data"][s] = dict(
            hgvs_errors["data"][s], **update
        )
    else:
        hgvs_errors["data"][s] = update

with open("hgvs_errors_new.json", "w") as hgvs_file:
    json.dump(hgvs_errors, hgvs_file, indent=4)
