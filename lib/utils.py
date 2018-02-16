'''
Utility functions for JSON Parsing
and
some pandas dataframe utilities implementing more exotic functions
'''

import re
import pandas

RE_HGVS = re.compile('[gcmnrpGCMNRP]\.')


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
