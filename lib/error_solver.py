import os
import csv
import json

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
    def __init__(self, hgvs_error_file='hgvs_errors.json', save=True):
        '''
        Params:
            hgvs_error_file Location to load the error-checking dict from and to.
            save Whether additions to the hgvs_errors should be automatically saved. This option is mainly useful for debugging purposes.
        '''
        self._hgvs_filepath = hgvs_error_file
        self._load()
        self._save_enabled = save

    def filename(self):
        return self._hgvs_filepath

    def _load(self):
        if os.path.exists(self._hgvs_filepath):
            j = json.load(open(self._hgvs_filepath, 'r'))
        else:
            j = {}
        self._data = j

    def __getitem__(self, entry):
        '''Genomic entry id as key
        '''
        if not entry:
            return entry
        key = str(entry['entry_id'])
        if key in self._data:
            return self._data[key]['cleaned']
        else:
            return entry

    def add_faulty(self, entry, wrong_strings):
        if not entry:
            return False
        key = str(entry['entry_id'])
        if key in self._data:
            return True
        else:
            self._data[key] = {'info' : entry,
                    'wrong' : wrong_strings,
                    'cleaned' : [] }
            if self._save_enabled:
                self._save()

    def _save(self):
        with open(self._hgvs_filepath, 'w') as f:
            json.dump(self._data, f, indent=4)
