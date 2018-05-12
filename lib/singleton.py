'''
Enable argumentless instantation and lazy configuration of API objects.
'''
from functools import wraps
from types import FunctionType

def force_configuration(func):
    @wraps(func)
    def force(*args, **kwargs):
        self = args[0]
        if not self.configured:
            raise RuntimeError("{} not configured yet.".format(self.__name__))
        return func(*args, **kwargs)
    return force


def configured(func):
    @wraps(func)
    def set_configure(*args, **kwargs):
        self = args[0]
        self.configured = True
        return func(*args, **kwargs)
    return set_configure


def init_conf(func):
    @wraps(func)
    def create_var(*args, **kwargs):
        self = args[0]
        self.configured = False
        return func(*args, **kwargs)
    return init_conf


class LazyConfigure:

    def __init__(self):
        self.configured = False

    def configure(self):
        self.configured = True

    def __new__(cls, *args, **kwargs):
        instance = super().__new__(cls, *args, **kwargs)
        inst_dict = getattr(instance, "__dict__")
        for attr_name, attr in inst_dict.items():
            if (type(attr) == "function" and (not attr_name.startswith("_")
                                              or attr_name == "configure")):
                inst_dict[attr_name] = force_configuration(attr)
        setattr(instance, "__dict__", inst_dict)

        return instance
