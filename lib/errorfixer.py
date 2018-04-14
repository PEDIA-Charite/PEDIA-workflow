'''
Error fixing in JSON files.
---
This class is mainly used to handle erroneous hgvs codes.
'''
import os
import json
from typing import Union, Tuple
from lib.constants import HGVS_ERRORDICT_VERSION

# Typing definitions for clearer type hints
ENTRY_ID = Union[str, int]
ERROR_ENTRY = Tuple[list, list, list]


class ErrorFixer:
    '''Map erroneous variant information to correct hgvs strings.
    '''
    def __init__(self,
                 hgvs_error_file: str = 'hgvs_errors.json',
                 hgvs_new_errors: str = 'hgvs_new_errors.json',
                 config: Union[None, "ConfigManager"] = None,
                 version: Union[None, int] = None):
        '''
        Params:
            hgvs_error_file: Location to load the error-checking info from.
            hgvs_new_errors: Location to send new errors to.
            save: Whether additions to the hgvs_errors should be automatically
                  saved. This option is mainly useful for debugging purposes.
            config: config object to directly load options
            version: override library internal version check
        '''
        if config:
            self._error_path = config.errorfixer["error_path"]
            self._new_error_path = config.errorfixer["new_error_path"]
        else:
            self._error_path = hgvs_error_file
            self._new_error_path = hgvs_new_errors

        self._error = self.load(self._error_path)
        latest_error_version = HGVS_ERRORDICT_VERSION \
            if version is None else version
        if self._error['version'] < latest_error_version:
            raise TypeError(
                "HGVS errordict version {} is older than {}.".format(
                    self._error["version"],
                    latest_error_version
                )
            )

        # always only capture errors in the newest iteration
        self._new = {}

    def get_filepath(self):
        '''Get path of the error file.'''
        return self._error_path

    def get_new_filepath(self):
        '''Get filepath of the new error file.'''
        return self._new_error_path

    def new_error(self, key: ENTRY_ID, value: ERROR_ENTRY) -> None:
        '''Add key and value to new error dictionary. These should later
        be manually checked and transferred to the correct hgvs errors dict.
        '''
        if key in self._new:
            self._new[key]['info'] += value[0]
            self._new[key]['correct'] += value[1]
            self._new[key]['wrong'] += value[2]
        else:
            self._new[key] = {
                'info': value[0],
                'correct': value[1],
                'wrong': value[2],
                'cleaned': []
            }
        self.save(self._new, self._new_error_path)

    def __getitem__(self, key: ENTRY_ID) -> dict:
        key = str(key)
        return self._error['data'][key]['cleaned']

    def get_data(self, key: ENTRY_ID):
        key = str(key)
        return self._error['data'][key]

    def __contains__(self, key: ENTRY_ID) -> bool:
        return str(key) in self._error['data']

    def __setitem__(self, key: ENTRY_ID, value: ERROR_ENTRY) -> bool:
        '''Add a faulty genomic entry and optionally wrong automatically
        generated hgvs strings. Both can be used as reference for the
        manual correction of hgvs strings.
        '''
        if not key:
            return False

        if not self.__contains__(key):
            self.new_error(key, value)

        return True

    @staticmethod
    def load(path: str) -> dict:
        '''Load an existing error json or initiate an empty dictionary.
        '''
        if os.path.exists(path):
            with open(path, "r") as hgvs_file:
                j = json.load(hgvs_file)
            # explicitly fail if error dict has no version
            if 'version' not in j:
                raise TypeError("HGVS error dict has no version field.")
        else:
            j = {"version": 0, "data": {}}
        return j

    @staticmethod
    def save(data, path):
        '''Save json data to a file with visual indentation.'''
        # do not save if path is ""
        if path:
            with open(path, 'w') as saved_dict:
                json.dump(data, saved_dict, indent=4)
