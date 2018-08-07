# Standard library
import os
import pprint
import timeit
import collections
import sys
# Third party library
import pandas as pd
# Local script
import alignment.commons as paths

__author__ = 'peeyush'

Path = paths.path()
basepath = Path.basepath


def filterpeaks(peak_data, name, filtering=True):
    # filtered_peak_data = {}
    # for k, v in peak_data.iteritems():
    # print k, v.shape
    name = name
    df = peak_data
    df['chr'] = df['chr'].astype('str')
    sample = name.split(' ')
    colnames = df.columns.values.tolist()
    if filtering:
        #print(colnames)
        indices1 = [i for i, s in enumerate(colnames) if sample[0] in s]
        indices2 = [i for i, s in enumerate(colnames) if sample[2] in s]
        # indices3 = [i for i, s in enumerate(colnames) if "Input" in s]
        for i in indices1:
            if ("RA" not in colnames[i]) and ("norm" not in colnames[i]) and ('Tag count' in colnames[i]):
                condition = colnames[i]
                break
            elif ("RA" in colnames[i]) and ("norm" not in colnames[i]) and ('Tag count' in colnames[i]):
                condition = colnames[i]
                break
        for i in indices2:
            if ("RA" not in colnames[i]) and ("norm" not in colnames[i]) and ('Tag count' in colnames[i]):
                control = colnames[i]
                break
            elif ("RA" in colnames[i]) and ("norm" not in colnames[i]) and ('Tag count' in colnames[i]):
                control = colnames[i]
                break
        else:
            raise ValueError("Filtering Sample name differs from column name.")
        print('Sample lane:' + condition)
        print('Control lane:' + control)
        # inputcol = colnames[indices3[0]]
        # print inputcol
        exclude_from_filtering = ['H3K36me3', 'H3K27me3', 'H3K4me1']

        ## condition for simple filtering of preaks
        if any(s in condition for s in exclude_from_filtering):
            print('Sample name in simple filtering list')
            df1 = df[df[condition] >= 2 * df[control]]
            final = df1
        else:
            print('Using default filtering....')
            df1 = df[df[condition] >= 2 * df[control]]
            df2 = df1[((df1['stop'] - df1['start']) / df1[condition]) <= 15]
            final = df2
        print('Default peak count:', df.shape)
        print('Filtered peak count:', final.shape)
    else:
        final = df
        print('Dataframe is not filtered:', final.shape)
    with open(basepath + '/further_analysis/filtered/filteredPeaksCount.txt', 'a') as file:
        file.write(name + '\t' + str(len(df)) + '\t' + str(len(final)) + '\n')
    # filtered_peak_data[name] = final
    dir_path = basepath + '/further_analysis/filtered/' + name
    paths.ensure_path(dir_path)
    samPath = os.path.join(dir_path, name + '.tsv')
    final.to_csv(samPath, sep="\t", header=True)
    final.index = range(len(final))
    return final, dir_path


def peaks_in_allsamples(peaksdirpath=None):
    '''
    Create a txt file for peaks count in all the files in csv folder
    :param peaksdirpath:
    :return:
    '''
    import pandas as pd
    if peaksdirpath is None:
        peaksdirpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/'
    # print(os.listdir(peaksdirpath))
    files = os.listdir(peaksdirpath)
    file = open(basepath + '/further_analysis/RawPeakCount.txt', 'w')
    file.write('Sample_name\tcontrol_sample\traw_peaks\n')
    if len(files) > 0:
        for peakfile in files:
            sample = peakfile.split('vs')
            if (os.path.isfile(os.path.join(peaksdirpath, peakfile))) and (peakfile.endswith('.csv')):
                # print('here')
                df = pd.read_csv(os.path.join(peaksdirpath, peakfile), sep='\t', header=0)
                if len(sample) > 1:
                    file.write(sample[0].strip() + '\t' + sample[1].strip() + '\t' + str(len(df)) + '\n')
                else:
                    file.write(sample[0].strip() + '\t--\t' + str(len(df)) + '\n')
    print('Done!!')
    file.close()


def extract_data(dataset, list_of_values, column):
    '''
    This will filter data based on values in column provided
    :param dataset: dataframe
    :param idlist: list
    :return: dataframe
    '''
    newdataset = dataset[dataset[column].isin(list_of_values)]
    return newdataset


def get_union_of_peaks(dict_of_peaks, outpath, filename=None):
    '''
    This function will calculate a union of peaks from given peaks dataframes.
    input = {'dataframe1': df1, 'datafram2':df2}
    :return: pandas.dataframe
    '''
    import os
    import sys
    import pprint
    import timeit
    import pandas as pd
    import collections
    df_dict = dict_of_peaks
    df_gr_dict = {}

    if filename is None:
        filename = 'Union_of_peaks_'+'_'.join(list(df_dict.keys()))

    # columns used for extracting data
    columns = ['chr', 'start', 'stop', 'length', 'summit', 'Next transcript strand',
               'GenomicPosition TSS=1250 bp, upstream=5000 bp',
               'Next Transcript stable_id', 'Next Transcript tss distance', 'Next transcript gene name',
               ]
    # Checking for required columns
    try:
        for name, df in df_dict.items():
            sys.stdout.write("Peaks in %s: %i \n" % (name, len(df)))
            df_dict[name] = df[columns]
    except Exception as e:
        print('Probably data frame does not have columns needed check below')
        pprint.pprint(columns)
        raise AttributeError(e)

    # Creating chromosome base grouping of peak df
    for name, df in df_dict.items():
        df_dict[name] = df[columns]
        df['chr'] = df['chr'].astype(str)
        df_gr_dict[name] = df.groupby('chr', sort=True)

    # finding overlap between dataframe
    start = timeit.default_timer()
    ignore_index = collections.defaultdict(list)  # creating dict to store overlap index for every df
    df_name = list(df_dict.keys())
    ignore_index[df_name[0]] = []  # Initialize first data frame as empty because all of the rows will be included

    for i_df in range(0, len(df_name)):
        peak_df_gr = df_gr_dict[df_name[i_df]]

        for i_df1 in range(i_df+1, len(df_name)):
            print('Comparing:', df_name[i_df], 'with', df_name[i_df1])
            peak_df1_gr = df_gr_dict[df_name[i_df1]]

            for gr, df in peak_df_gr:
                if gr in peak_df1_gr.groups:
                    df1 = peak_df1_gr.get_group(gr)
                    for ind, row in df.iterrows():
                        for ind1, row1 in df1.iterrows():
                            if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                                ignore_index[df_name[i_df1]].append(ind1)
    #print(ignore_index.keys())

    # last step joining data frames after removing overlap peaks
    unique_df = pd.DataFrame(columns=columns)

    for key, index2ignore in ignore_index.items():
        #print(df)
        df = df_dict[key][columns]
        #print(key, set(index2ignore), len(set(index2ignore)))
        df = df[~df.index.isin(set(index2ignore))]
        df['parent_df'] = key
        unique_df = unique_df.append(df, ignore_index=True)

    columns.append('parent_df')
    unique_df = unique_df.reindex(columns=columns)
    print('Combined peaks count:', len(unique_df))
    unique_df.to_csv(os.path.join(outpath, filename + '.tsv'), header=True, sep='\t', index=True)
    stop = timeit.default_timer()
    print('Time consumed in analysis:', stop-start/60)
    return unique_df
