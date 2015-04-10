__author__ = 'peeyush'

import pandas as pd

def extract_keys(path, dataframe, columnNames):
    File = open(path)
    keys = File.readline()
    File.close()
    return filter_rows(keys, dataframe, columnNames)

def filter_rows(keys, dataFrame, columnName):
    '''
    This Method will filter a dataframe using list of keys and column names
    :param keys:
    :param dataFrame:
    :param columnName:
    :return:
    '''
    new_df = pd.DataFrame()
    for key in keys:
        for index, rows in dataFrame.iterrows():
            if rows[columnName] == key:
                new_df = new_df.append(rows)
    return new_df

