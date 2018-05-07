'''
Functions to handle quality check procedures.

Quality check logs follow identical structure.

section - {case_id: data}
'''
import json
import os


def read_json(path: str):
    '''Read a json file.'''
    with open(path) as jsfile:
        data = json.load(jsfile)
    return data


def read_old_new(path: str):
    '''Read old and new file.
    Old files are suffixed with .old. If no old file can be found the current
    new file will be returned.
    '''

    new = read_json(path)

    oldpath = path+".old"
    if os.path.exists(oldpath):
        old = read_json(oldpath)
    else:
        old = new

    return new, old


def diff_keys(newdata: dict, olddata: dict) -> dict:
    '''Compare existence of keys in new and old log.'''
    newkeys = set(newdata.keys())
    oldkeys = set(olddata.keys())

    old_only = oldkeys - newkeys
    new_only = newkeys - oldkeys

    return {
        "old_only": list(old_only),
        "new_only": list(new_only)
    }


def diff_log(newlog: dict, oldlog: dict) -> dict:
    '''Compare results of two logs.
    '''

    diff_summary = {
        section: diff_keys(
            newdata, oldlog[section] if section in oldlog else {}
        )
        for section, newdata in newlog.items()
    }
    return diff_summary


def diff_quality_check(config_data: "ConfigParser"):
    '''Get difference of QC results in comparison to last run.'''
    # get json qc files
    json_qc_path = config_data.jsonparser["json_qc_log"]
    qc_path = config_data.quality["qc_detailed_log"]

    json_qc_diff = diff_log(*read_old_new(json_qc_path))
    qc_diff = diff_log(*read_old_new(qc_path))

    all_diff = {**json_qc_diff, **qc_diff}

    print("== Quality check difference to last run ==")

    print(
        "\n".join(
            [
                "= {} =\n{}".format(
                    section,
                    "\n".join([
                        "{}: {}".format(k, ", ".join(v))
                        for k, v in diff.items()
                    ])
                )
                for section, diff in all_diff.items()
            ]
        )
    )
