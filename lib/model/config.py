from configparser import ConfigParser

class Indexer:
    def __init__(self, obj, section):
        self._obj = obj
        self._section = section

    def __getitem__(self, key):
        return self._obj.get(self._section, key)

class ConfigManager:
    """Configuration for all steps in the pipeline.
    """
    def __init__(self, conf_file='config.ini'):
        self.load_config(conf_file)

        # adding each section as an attribute to the class
        for section in self._data.sections():
            setattr(self, section, Indexer(self, section))

    def load_config(self, conf_file):
        """Load necessary vars from config file.
        """
        self._data = ConfigParser()
        self._data.read(conf_file)

    def get(self, section, key):
        if isinstance(key, tuple) or isinstance(key, list):
            return {k:self.get(section, k) for k in key}
        elif isinstance(key, str):
            return self._data[section][key]
        else:
            raise TypeError

    def __getitem__(self, key):
        return self._data[key]
