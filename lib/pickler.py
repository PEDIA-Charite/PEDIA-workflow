import pickle

from lib.model.config import PEDIAConfig


class CasePickler(pickle.Pickler):

    def persistent_id(self, obj):
        if isinstance(obj, PEDIAConfig):
            return ("PEDIAConfig", obj.path)
        else:
            return None


class CaseUnpickler(pickle.Unpickler):

    def __init__(self, filedesc):
        super().__init__(filedesc)
        self.configs = {}

    def persistent_load(self, pid):
        type_tag, key_id = pid
        if type_tag == "PEDIAConfig":
            if key_id not in self.configs:
                self.configs[key_id] = PEDIAConfig(key_id)
            return self.configs[key_id]
        else:
            raise pickle.UnpicklingError("unsupported persistent object")
