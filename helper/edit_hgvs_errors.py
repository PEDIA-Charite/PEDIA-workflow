#!/bin/env python3
from argparse import ArgumentParser
import json
import os
import sys


def exit():
    print("Exiting...")
    sys.exit(0)


def print_info(errordict: dict):
    '''Print some information on current hgvs dict.'''
    print("HGVS errordict version", errordict["version"])
    print(len(errordict["data"]), "entries")


def parse_args():
    parser = ArgumentParser(
        description=(
            "Enter new entries into the hgvs errors.json, "
            "while checking the current entry."
        )
    )
    parser.add_argument("hgvs_file")
    args = parser.parse_args()
    return args


def tokenize(cmd):
    return cmd.split()


class HGVSEditor:
    def __init__(self, args):
        with open(args.hgvs_file, "r") as hgfile:
            oldhgvs = json.load(hgfile)
        self.old = oldhgvs
        self.newdata = {}

    def run_browser(self):
        print_info(self.old)
        while True:
            self.prompt()

    def prompt(self):
        try:
            cmd = input("> ")
        except EOFError:
            exit()
        except KeyboardInterrupt:
            exit()

        toks = tokenize(cmd)
        if toks[0] == "help":
            print("Usage: ls, show, new")
        elif toks[0] == "ls":
            print("\n".join(self.old["data"].keys()))
        elif toks[0] == "show":
            if len(toks) != 2:
                print("Usage: show <entry_id>")
            print(json.dumps(self.by_id(toks[1]), indent=4))
        elif toks[0] == "new":
            if len(toks) != 2:
                print("Usage: new <entry_id>")
            self.new_entry(toks[1])
        elif toks[0] == "edit":
            if len(toks) != 2:
                print("Usage: edit <entry_id>")
            self.edit_entry(toks[1])
        elif toks[0] == "save":
            self.save_entries
        elif toks[0] == "changes":
            print(json.dumps(self.newdata, indent=4))

    def by_id(self, search):
        if search in self.old["data"]:
            return self.old["data"][search]
        elif search in self.newdata:
            return self.newdata[search]
        else:
            return {"error": "{} not found".format(search)}

    def new_entry(self, entry_id):
        if entry_id in self.old:
            print("Entry already exists. Update or delete first.")
            return
        print("CREATING new entry")
        note = input("Enter note: ")
        cleaned = [
            s.strip() for s in input("Enter cleaned (sep ;) :").split(";")
        ]
        self.newdata[entry_id] = {
            "note": note,
            "cleaned": cleaned,
        }

    def edit_entry(self, entry_id):
        old = self.by_id(entry_id)
        if "error" in old:
            print(old["error"])
        print("EDITING", entry_id)
        print(json.dumps(old, indent=4))
        print("Old note: ", old["note"] if "note"in old else "")
        note = input("Enter new(enter to skip)> ")
        if note:
            old["note"] = note
        print("Old cleaned: ", old["cleaned"] if "cleaned"in old else "")
        cleaned = input("Enter new cleaned(enter to skip)> ")
        if cleaned:
            old["cleaned"] = cleaned

        print("Saving edits.")
        self.newdata[entry_id] = old

    def save_entries(self):
        dumpdata = {
            "version": self.old["version"] + 1,
            "data": {**self.old["data"], **self.newdata}
        }

        new_path = "hgvs_errors_edited.json"
        with open(new_path, "r") as hgnew:
            json.dump(dumpdata, hgnew, indent=4)

        print("HGVS error dict saved to {}".format(new_path))


def main():
    editor = HGVSEditor(parse_args())
    editor.run_browser()


if __name__ == "__main__":
    main()
