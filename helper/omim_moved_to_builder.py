import re
import json
import pandas

mim_titles = pandas.read_table("data/mimTitles.txt", sep="\t", comment="#",
                               names=["prefix", "mim_number", "title", "alt",
                                       "inc"])
removed = mim_titles.loc[mim_titles["title"].str.contains("REMOVED FROM DATABASE")]

moved_to = mim_titles.loc[mim_titles["title"].str.contains("MOVED TO")]
moved_dict = dict(zip(moved_to["mim_number"], moved_to["title"]))

omim_trans = {r: [] for r in removed["mim_number"]}

for mim, msg in moved_dict.items():
    matches = re.findall("\d+", msg)
    omims = [m for m in matches]
    if mim in omim_trans:
        raise RuntimeError("Lol no")
    omim_trans[mim] = omims

with open("omim_deprecated_replacement.txt", "w") as repfile:
    json.dump(omim_trans, repfile, indent=4)
