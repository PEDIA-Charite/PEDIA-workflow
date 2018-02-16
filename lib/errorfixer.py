'''
Error fixing in JSON files.
---
This class is mainly used to handle erroneous hgvs codes.
'''
import os
import json
from typing import Union


class ErrorFixer:
    '''Map erroneous variant information to correct hgvs strings.
    These functions are a last resort for manual intervention.
    The interventions should be used like this:
    genomic_entries ids are mapped to a list of hgvs codes
    the errordict has the following structure:
    {entry_id : {
        info:<current info as reference here>,
        cleaned:[<hgvs strings>]
        }}
    '''
    def __init__(self,
                 hgvs_error_file: str='hgvs_errors.json',
                 save: bool=True):
        '''
        Params:
            hgvs_error_file: Location to load the error-checking dict from and
                             to.
            save: Whether additions to the hgvs_errors should be automatically
                  saved. This option is mainly useful for debugging purposes.
        '''
        self._hgvs_filepath = hgvs_error_file
        self._load()
        self._save_enabled = save

    def get_filepath(self):
        return self._hgvs_filepath

    def _load(self):
        '''Load an existing error json or initiate an empty dictionary.
        '''
        if os.path.exists(self._hgvs_filepath):
            j = json.load(open(self._hgvs_filepath, 'r'))
        else:
            j = {}
        self._data = j

    def __getitem__(self, entry: Union[str, int]) -> [str]:
        '''Genomic entry id as key and return a list of cleaned hgvs codes.
        '''
        if not entry:
            return entry
        key = str(entry['entry_id'])
        if key in self._data:
            return self._data[key]['cleaned']
        else:
            return entry

    def add_faulty(self, entry: dict, wrong_strings: [str]) -> bool:
        '''Add a faulty genomic entry and optionally wrong automatically
        generated hgvs strings. Both can be used as reference for the
        manual correction of hgvs strings.
        '''
        if not entry:
            return False
        key = str(entry['entry_id'])
        if key in self._data:
            return True
        else:
            self._data[key] = {'info': entry,
                               'wrong': wrong_strings,
                               'cleaned': []}
            if self._save_enabled:
                self._save()

    def _save(self):
        with open(self._hgvs_filepath, 'w') as saved_error_dict:
            json.dump(self._data, saved_error_dict, indent=4)
