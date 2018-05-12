'''
Functions to handle quality check procedures.

Quality check logs follow identical structure.

section - {case_id: data}
'''
import json
import os
import shutil

from lib.visual import progress_bar


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


def diff_quality_check(qc_path):
    '''Get difference of QC results in comparison to last run.'''
    # get json qc files
    all_diff = diff_log(*read_old_new(qc_path))

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


class QualityCheck:

    def __init__(self, args, processor):

        self._qc_logs = None

        self.data = self.processor

        self._write_logs = not args.single

        self.qc_logpath = self.data.config.quality["qc_detailed_log"]

    def summary(self):
        qc_passed = self.qc_logs["passed"]
        qc_failed_msg = self.qc_logs["failed"]
        qc_vcf_failed = self.qc_logs["vcf_failed"]
        return {
            "pass": len(qc_passed),
            "fail": len(qc_failed_msg) + len(qc_vcf_failed)
        }, qc_passed

    def _create_logs(self):
        all_logs = {**self.qc_jsons(), **self.qc_cases()}
        return all_logs

    @property
    def qc_logs(self):
        if not self._qc_logs:
            self._qc_logs = self._create_logs()
        return self._qc_logs

    def _write_log(self):
        logpath = self.qc_logpath

        if os.path.exists(logpath):
            shutil.move(logpath, logpath+".old")

        with open(logpath, "w") as failedfile:
            json.dump(self.qc_logs, failedfile, indent=4)


    def qc_jsons(self):
        qc_failed_results = {
            "json_check_failed": {
                case_id: {
                    "issues": issues,
                }
                for case_id, (valid, issues) in
                [(j.get_case_id(), j.check()) for j in self.data.jsons]
                if not valid
            }
        }
        return qc_failed_results

    def qc_cases(self):
        '''Output quality check summaries.'''

        qc_cases = self.data.qc_cases

        # Cases failing qc altogether
        qc_failed_msg = {
            case.case_id: valid for valid, case in qc_cases if not valid[0]
        }
        # cases with multiple diagnosis
        multi_no_omim = {
            k: v for k, v in
            [(c.case_id, c.get_diagnosis()) for _, c in qc_cases]
            if any(str(d["omim_id"]) == "0" for d in v) and len(v) > 1
        }

        # Cases passing quality check
        qc_passed = {
            case.case_id: case for valid, case in qc_cases if valid[0]
        }

        qc_vcf = {
            case_id: case.check_vcf() for case_id, case in qc_passed.items()
        }

        qc_vcf_failed = {
            case_id: result
            for case_id, result in qc_vcf.items() if not result[0]
        }

        # cases have to pass vcf check
        qc_passed = {
            case_id: case
            for case_id, case in qc_passed.items() if qc_vcf[case_id][0]
        }

        # Cases with mutations marked as benign excluded from analysis
        qc_benign_passed = {
            k: v for k, v in
            {k: v.get_benign_excluded() for k, v in qc_passed.items()}.items()
            if v > 0
        }

        # Cases where pathogenic diagnosed mutation is not in geneList
        @progress_bar("Get pathogenic genes in geneList")
        def pathogenic_genes_process(cases):
            '''Get boolean value, whether pathogenic gene is contained in
            genes converted from detected syndromes.'''
            for case_id, case_obj in cases.items():
                yield case_id, case_obj.pathogenic_gene_in_gene_list(
                    self.data.config.omim
                )

        qc_pathongenic_passed = dict(
            c for c in pathogenic_genes_process(qc_passed) if not c[1][0]
        )

        # Compiled stats to be dumped into a json file
        qc_output = {
            "failed": qc_failed_msg,
            "benign_excluded": qc_benign_passed,
            "pathogenic_missing": qc_pathongenic_passed,
            "vcf_failed": qc_vcf_failed,
            "multi_no_omim": multi_no_omim,
            "passed": {k: '' for k in qc_passed.keys()},
        }

        if not self._write_logs:
            print(json.dumps(qc_output, indent=4))

        return qc_output
