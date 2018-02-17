'''
Error fixing in JSON files.
---
This class is mainly used to handle erroneous hgvs codes.
'''
import os
import json
from typing import Union, Callable

# Typing definitions for clearer type hints
ENTRY_ID = Union[str, int]


def save_content(func: Callable):
    '''Decorator for saving the dictionary.
    '''
    def wrapper(*args, **kwargs):
        '''Save the dictionary after executing the function.
        '''
        self = args[0]
        resp = func(*args, **kwargs)
        self._parent._save(self._obj, self._path)
        return resp
    return wrapper


class Index:
    '''Base index class.
    '''
    def __init__(self, parent: 'ErrorFixer', obj: dict, path: str):
        self._obj = obj
        self._parent = parent
        self._path = path

    def __getitem__(self, key: ENTRY_ID) -> dict:
        key = str(key)
        return self._obj[key]

    def __contains__(self, key: ENTRY_ID) -> bool:
        return str(key) in self._obj

    @save_content
    def __setitem__(self, key: ENTRY_ID, value: dict):
        self._obj[key] = value


class ErrorIndex(Index):
    '''Error index class.
    '''
    def __getitem__(self, key: ENTRY_ID) -> list:
        '''Genomic entry id as key and return a list of cleaned hgvs codes.
        '''
        error_dict = super().__getitem__(key)
        return error_dict['cleaned']

    def __setitem__(self, key: ENTRY_ID, value: (dict, list, list)) -> bool:
        '''Add a faulty genomic entry and optionally wrong automatically
        generated hgvs strings. Both can be used as reference for the
        manual correction of hgvs strings.
        '''
        if not key:
            return False
        if self.__contains__(key):
            return True
        else:
            value = {'info': value[0],
                     'wrong': value[1],
                     'cleaned': value[2]
                     }
            super().__setitem__(key, value)
            return True


class PartialIndex(Index):
    '''Partial index class.
    '''
    def __getitem__(self, key):
        error_dict = super().__getitem__(key)
        return error_dict

    def __setitem__(self, key: ENTRY_ID, value: (list, list)) -> bool:
        '''Add partially faulty hgvs strings for manual confirmation, that
        no relevant mutations have been missed.
        Manually correcting these, will require adding an entry to the error
        dictionary.
        '''
        value = {'correct': value[0], 'wrong': value[1]}
        super().__setitem__(key, value)


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
    Files with multiple hgvs candidates, where errors are generated will be
    saved to the partial dictionary, which is only used for verification
    purposes.
    If keys have already been added to the partial dictionary, messages about
    hgvs parsing errors will be suppressed.
    '''
    def __init__(self,
                 hgvs_error_file: str='hgvs_errors.json',
                 hgvs_partial_file: str='hgvs_partial.json',
                 save: bool=True):
        '''
        Params:
            hgvs_error_file: Location to load the error-checking dict from and
                             to.
            save: Whether additions to the hgvs_errors should be automatically
                  saved. This option is mainly useful for debugging purposes.
        '''
        self._error_path = hgvs_error_file
        self._partial_path = hgvs_partial_file
        self._load('_error', self._error_path)
        self._load('_partial', self._partial_path)
        self._save_enabled = save

        self.error = ErrorIndex(self, self._error, self._error_path)
        self.partial = PartialIndex(self, self._partial, self._partial_path)

    def get_filepath(self):
        return self._error_path

    def get_partial_filepath(self):
        return self._partial_path

    def _load(self, attrname, path):
        '''Load an existing error json or initiate an empty dictionary.
        '''
        if os.path.exists(path):
            j = json.load(open(path, 'r'))
        else:
            j = {}
        setattr(self, attrname, j)

    def _save(self, data, path):
        with open(path, 'w') as saved_dict:
            json.dump(data, saved_dict, indent=4)
