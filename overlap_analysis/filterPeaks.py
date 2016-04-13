__author__ = 'peeyush'
import alignment.commons as paths
import os
Path = paths.path()
basepath = Path.basepath


def filterpeaks(peak_data, name, filtering=True):
    #filtered_peak_data = {}
    file = open(basepath + '/further_analysis/filtered/filteredPeaksCount.txt', 'w')
    file.write('Sample_name\traw_peaks\tfiltered_peaks\n')
    #for k, v in peak_data.iteritems():
    #print k, v.shape
    name = name
    df = peak_data
    sample = name.split(' ')
    colnames = df.columns.values.tolist()
    indices1 = [i for i, s in enumerate(colnames) if sample[0] in s]
    indices2 = [i for i, s in enumerate(colnames) if sample[2] in s]
    #indices3 = [i for i, s in enumerate(colnames) if "Input" in s]
    for i in indices1:
        if "RA" not in colnames[i] and "norm" not in colnames[i]:
            condition = colnames[i]
            break
        elif "RA" in colnames[i] and "norm" not in colnames[i]:
            condition = colnames[i]
            break
    for i in indices2:
        if "RA" not in colnames[i] and "norm" not in colnames[i]:
            control = colnames[i]
            break
        elif "RA" in colnames[i] and "norm" not in colnames[i]:
            control = colnames[i]
            break
    else:
        raise ValueError("Filtering Sample name differs from column name.")
    print control
    print condition
    #inputcol = colnames[indices3[0]]
    #print inputcol
    if filtering:
        ## condition for simple filtering of preaks
        df1 = df[df[condition] >= 3*df[control]]
        df2 = df1[((df1['stop']-df1['start'])/df1[condition]) <= 15]
        final = df2
        print final.shape
    else:
        final = df
        print final.shape
    file.write(name.split(' ')[0]+'\t'+str(len(df))+'\t'+str(len(final))+'\n')
    #filtered_peak_data[name] = final
    dirPATH = basepath + '/further_analysis/filtered/'+name
    paths.ensure_path(dirPATH)
    samPath = os.path.join(dirPATH, name+'.txt')
    final.to_csv(samPath, sep="\t", header=True)
    file.close()
    return final, dirPATH


def extract_data(dataset, list_of_values, column):
    '''
    This will filter data based on values in column provided
    :param dataset: dataframe
    :param idlist: list
    :return: dataframe
    '''
    newdataset = dataset[dataset[column].isin(list_of_values)]
    return newdataset