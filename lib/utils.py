import re

def optional_descent(data, key_trail, default=''):
    '''Follow a iterable object of key names in a nested dict.'''
    d = data
    # do not attempt if not a dict
    if data == '':
        return data

    for k in key_trail:
        if k in d:
            d = d[k]
        else:
            return default
    return d

RE_HGVS = re.compile('[gcmnrpGCMNRP]\.')
def checkHGVS(string):
    # trick to remove all whitespace
    string = "".join(string.split())
    m = RE_HGVS.search(string)
    return m is not None

def list_all_in(a, b):
    '''Check that all elements of a are in b.'''
    if len(a) != len(b):
        return False
    for i in a:
        if i not in b:
            return False
    return True
