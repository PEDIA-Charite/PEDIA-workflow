'''
Utility functions for JSON Parsing
and
some pandas dataframe utilities implementing more exotic functions
'''
import re
import hashlib
from typing import Union

import pandas


RE_HGVS = re.compile(r'[gcmnrpGCMNRP]\.')


def optional_descent(data: dict, key_trail: [str], default: str = '') \
        -> Union[dict, list, str]:
    '''Follow a iterable object of key names in a nested dict.'''
    if data == '':
        return data

    for k in key_trail:
        if k in data:
            data = data[k]
        else:
            return default
    return data


def check_hgvs(hgvs_string: str):
    '''Check whether a possible HGVS String might be contained inside a string
    using RE_HGVS, which only checks expressions of the form 'x.<etc>'
    '''
    # remove all whitespace
    string = "".join(hgvs_string.split())
    match = RE_HGVS.search(string)
    return match is not None


def list_all_in(tested: list, reference: list) -> bool:
    '''Check that all elements of tested list are in reference list.'''
    if len(tested) > len(reference):
        return False
    for i in tested:
        if i not in reference:
            return False
    return True


def explode_df_column(dataframe: pandas.DataFrame, column: str) \
        -> pandas.DataFrame:
    '''Create multiple rows from a list inside a column.
    Args:
        df: Pandas dataframe to be exploded. This process should not change the
            object in place.
        column: Column loc identifier containing the list, which can be
                exploded.

    Returns:
        Dataframe with the selected list column as multiple rows. All other
        values will be duplicated.
    '''

    exploded = dataframe.apply(
        lambda x: pandas.Series(x[column]), axis=1).stack().reset_index(
            level=1, drop=True)
    exploded.name = column
    dataframe = dataframe.drop(
        [column], axis=1).join(exploded).dropna(subset=[column])
    return dataframe


def get_file_hash(filepath: str) -> str:
    '''Get MD5 Hash of file.'''
    md5hash = hashlib.md5()
    with open(filepath, "rb") as fileobj:
        while True:
            bytestr = fileobj.read(128)
            if not bytestr:
                break
            md5hash.update(bytestr)
    md5 = md5hash.hexdigest()
    return md5
