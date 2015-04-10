__author__ = 'peeyush'

from pandas import read_csv
import pandas as pd


df1 = read_csv('/home/peeyush/Desktop/HW_diff_exp/Exp1.csv', sep=',')
df2 = read_csv('/home/peeyush/Desktop/HW_diff_exp/Exp2.csv', sep=',')
df3 = read_csv('/home/peeyush/Desktop/HW_diff_exp/Exp3.csv', sep=',')
column = ['Affy_ID', 'FC_48', 'Diff_call_48', 'FC_72', 'Diff_call_72', 'FC_14d', 'Diff_call_14d']
new_df = pd.DataFrame()

for index, row in df1.iterrows():
    if '~' in row[8]:
        row[8] = row[8][1:]
    if '~' in row[10]:
        row[10] = row[10][1:]
    if '~' in row[16]:
        row[16] = row[16][1:]
for index, row in df2.iterrows():
    if '~' in row[8]:
        row[8] = row[8][1:]
    if '~' in row[10]:
        row[10] = row[10][1:]
    if '~' in row[16]:
        row[16] = row[16][1:]
for index, row in df3.iterrows():
    if '~' in row[8]:
        row[8] = row[8][1:]
    if '~' in row[10]:
        row[10] = row[10][1:]
    if '~' in row[16]:
        row[16] = row[16][1:]

for i in range(0, len(df1)-1):
    if df1.iloc[i][0] == df2.iloc[i][0] and df2.iloc[i][0] == df3.iloc[i][0]:
        Id = df1.iloc[i][0]

        FC_48 = (float(df1.iloc[i][8]) + float(df2.iloc[i][8]) + float(df3.iloc[i][8]))/3
        Diff_call_48 = df1.iloc[i][7] + df2.iloc[i][7] + df3.iloc[i][7]

        FC_72 = (float(df1.iloc[i][10]) + float(df2.iloc[i][10]) + float(df3.iloc[i][10]))/3
        Diff_call_72 = df1.iloc[i][9] + df2.iloc[i][9] + df3.iloc[i][9]

        FC_14d = (float(df1.iloc[i][16]) + float(df2.iloc[i][16]) + float(df3.iloc[i][16]))/3
        Diff_call_14d = df1.iloc[i][15] + df2.iloc[i][15] + df3.iloc[i][15]
        new_df = new_df.append(pd.Series([Id, FC_48, Diff_call_48, FC_72, Diff_call_72, FC_14d, Diff_call_14d]), ignore_index=True)

new_df.columns = column
new_df.to_csv('/home/peeyush/Desktop/HW_diff_exp/FC_table.csv', sep=',', header=True, index_label=False, index=False)

