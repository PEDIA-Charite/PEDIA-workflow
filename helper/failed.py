import json
import os

with open("hgvs_errors_v7.json", "r") as hgfile:
    olddata = json.load(hgfile)

dat = olddata["data"]

failedjannovar_corrected = [k for k, v in dat.items() if "note" in v and v["note"] == "failedjannovar_corrected"]

for k in failedjannovar_corrected:
    with open("process/aws_dir/cases/{}.json".format(k), "r") as jsfile:
        jsdata = json.load(jsfile)
        ge = jsdata["genomic_entries"]
    if len(ge) != 1:
        print("Too many entries", k, ge)
        continue
    if ge[0] not in dat:
        dat[ge[0]] = dat.pop(k)
    else:
        del dat[k]


olddata["version"] = 8

with open("hgvs_errors_v8.json", "w") as hgfile:
    json.dump(olddata, hgfile, indent=4)
