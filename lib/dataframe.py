import pandas

'''
Some pandas dataframe utilities implementing more exotic functions
'''
def explode_df_column(df, column):
    c = df.apply(lambda x: pandas.Series(x[column]), axis=1).stack().reset_index(level=1, drop=True)
    c.name = column
    df = df.drop([column], axis=1).join(c).dropna(subset=[column])
    return df
