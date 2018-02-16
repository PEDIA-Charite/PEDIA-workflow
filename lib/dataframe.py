'''
Some pandas dataframe utilities implementing more exotic functions
'''
import pandas


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
