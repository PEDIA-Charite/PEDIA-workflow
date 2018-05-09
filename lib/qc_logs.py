'''Load QC logs and extract information.'''
import os
import json

from lib.utils import load_json


class QCLogs:
    '''Concatenate QC logs and expose functions for simple case id queries.
    '''

    def __init__(self, log_paths: [str]):
        self.data = {
            k: v
            for d in [load_json(p, {}) for p in log_paths]
            for k, v in d.items()
        }

    def get_case_info(self, case_id: str) -> [dict]:
        '''Search a case id and return mentions in all sections.'''

        mentions = {
            section: data[case_id]
            for section, data in self.data.items()
            if case_id in data
        }

        return mentions
