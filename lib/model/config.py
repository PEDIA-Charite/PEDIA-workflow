'''
Configuration loader to get parameters from a ini file.
---
This is mainly used to define download directories and API access.
'''
from typing import Union, Iterable, Dict
from configparser import ConfigParser


class Indexer:
    '''Implement config sections as class attributes.
    '''
    def __init__(self, obj: 'ConfigManager', section: str):
        self._obj = obj
        self._section = section

    def getboolean(self, key: str):
        return self._obj[self._section].getboolean(key)

    def __getitem__(self, key: Union[str, tuple, list]):
        return self._obj.get(self._section, key)


class ConfigManager:
    """Configuration for all steps in the pipeline.
    """
    def __init__(self, conf_file: str = 'config.ini'):
        self._load_config(conf_file)

        # adding each section as an attribute to the class
        for section in self._data.sections():
            setattr(self, section, Indexer(self, section))

    def _load_config(self, conf_file: str):
        """Load necessary vars from config file.
        """
        self._data = ConfigParser()
        self._data.read(conf_file)

    def get(self, section: str, key: Union[str, Iterable[str]]) \
            -> Union[Dict[str, str], str, int, bool]:
        '''Get a single parameter or a list of parameters by name.
        If an iterable is given a dictionary containing key value pairs will be
        returned.
        '''
        if isinstance(key, tuple) or isinstance(key, list):
            return {k: self.get(section, k) for k in key}
        elif isinstance(key, str):
            return self._data[section][key]
        else:
            raise TypeError

    def __getitem__(self, key):
        return self._data[key]

    def __contains__(self, key):
        return key in self._data
