__author__ = 'peeyush'

import pandas as pd


class Overlaps():
    """
    This class object will hold filteredPeaks data and peaksList for calculation differential binding.
    """

    def __init__(self, overlappinglist, filter_peaks):
        self.samples_names = overlappinglist
        self.filter_peaks = filter_peaks


def diffBinding(self, basepeakfile):
    '''
    This function will extract summit (+-500) peak data if peak length is >1000 from provided peaks.
    This dataframe can be used with DESeq for differential bound calculation.
    :return:
    '''
    import pysam
    import sys
    print "\nCheck point: diffBinding"
    sample_name = self.samples_names
    dataframes = self.filter_peaks
    # print type(sample_name)
    #df = dataframes.get(sample_name[0]).iloc[:, 0:1]
    df = pd.DataFrame()
    df = pd.concat([df, dataframes.get(basepeakfile)['chr']], axis=1)
    df = pd.concat([df, dataframes.get(basepeakfile)['start']], axis=1)
    df = pd.concat([df, dataframes.get(basepeakfile)['stop']], axis=1)
    df = pd.concat([df, dataframes.get(basepeakfile)['GenomicPosition TSS=1250 bp, upstream=5000 bp']], axis=1)
    df = pd.concat([df, dataframes.get(basepeakfile)['Next transcript gene name']], axis=1)
    df = pd.concat([df, dataframes.get(basepeakfile)['Next transcript strand']], axis=1)
    df = pd.concat([df, dataframes.get(basepeakfile)['summit']], axis=1)
    df['cookiecut_start'] = 0
    df['cookiecut_stop'] = 0
    #print df.head()
    print df.dtypes
    for sample in sample_name:
        df[sample] = 0
        #print '\n'+sample
        sample_bam_path = getBam(sample) #.split(' vs ')[0]
        sample_bam = pysam.Samfile(sample_bam_path, "rb")
        for k, v in df.iterrows():
            sys.stdout.write("\rNumber of peaks processed:%d" % k)
            sys.stdout.flush()
            #print v['start'], v['summit']
            if v['stop']-v['start'] > 1000:
                chr = str(v['chr'])
                summit = v['start']+v['summit']
                tags = sample_bam.count(chr, summit-500, summit+500)
                df.loc[k,'cookiecut_start'] = summit-500
                df.loc[k,'cookiecut_stop'] = summit+500
                df.loc[k,sample] = tags
            else:
                chr = str(v['chr'])
                tags = sample_bam.count(chr, v['start'], v['stop'])
                df.loc[k,'cookiecut_start'] = v['start']
                df.loc[k,'cookiecut_stop'] = v['stop']
                df.loc[k,sample] = tags
        sample_bam.close()
    df.to_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/' + '_'.join(sample_name) + '.csv',
            sep=",", encoding='utf-8', ignore_index=True)

def getBam(name):
    '''
    This function takes the sample name and return the path of its associated bam file.
    :param name:
    :return:
    '''
    from os import listdir
    path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/'
    bam_list = listdir('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane')
    Dir = None
    file = None
    for i in bam_list:
        if name in i and 'dedup' in i:
            if 'RA' in name and 'RA' in i:
                Dir = path+i
                #print Dir
                for j in listdir(Dir):
                    if j.endswith('.bam'):
                        file = j
                        print '\nBam file selected: '+j
            if 'RA' not in name and 'RA' not in i:
                Dir = path+i
                #print Dir
                for j in listdir(Dir):
                    if j.endswith('.bam'):
                        file = j
                        print '\nBam file selected: '+j
    if file is None:
        raise KeyError('Bam file cannot be found for '+name)
    else:
        return Dir+'/'+file




def factor_seperate(dataframe, factor):
    """
    This function will sort the dataframe and return reletive positions of factor (eg. chromosome, genomic region) in peak table.
    :param overlapsObj:
    :return:
    """
    #if type(factor) != str:
    #    raise ValueError('def chr_position wrong column name')
    print 'Defined column is: '+factor
    peak_list = dataframe.sort([factor], ascending=True)
    # print df.head()
    chr_pos_list = peak_list[factor]
    chr_pos_list = chr_pos_list.astype(basestring)
    start_pos = 0
    index_counter = 0
    chr_pos = {}
    last_entity = None
    for i in chr_pos_list:
        if index_counter == 0:
            last_entity = i
            print last_entity
        index_counter += 1
        if i != last_entity:
            #print start_pos
            #print last_entity
            chr_pos[last_entity] = start_pos, index_counter - 1
            start_pos = index_counter
            last_entity = i
            #print last_entity
    chr_pos[last_entity] = start_pos, index_counter
    print 'Factor position: ', chr_pos
    return divide_dataframe(peak_list, chr_pos)

def divide_dataframe(dataframe, posDict):
    '''
    This method takes result from factor_seperate and divide dataframe into respective parts
    :param dataframe:
    :param posDict:
    :return:
    '''
    dict_df = {}
    for k, v in posDict.iteritems():
        sub_df = dataframe[v[0]:v[1]]
        dict_df[k] = sub_df
    return dict_df


def col_position(dataframe, colname):
    """
    Returns the column indices of a column name in df
    :param dataframe:
    :param colname:
    :return:
    """
    colname = colname.split(' ')[0]
    collist = dataframe.columns.values.tolist()
    indices = [i for i, s in enumerate(collist) if colname in s]
    return indices


def non_overlapping_peaks(self, overlap_no):
    """
    Reterive non overlapping data from peak data
    :param overlapingDF:
    :param df1:
    :param df2:
    :return: tuple of non-overlapping data
    """
    print "\nCheck point: Differential peaks selection"
    iter_overlaps = self.db_overlaps
    dataframes = self.filter_peaks

    for name, overlap in iter_overlaps.iteritems():
        n1 = name.split("_vs_")[0]  # split to get the sample name
        n2 = name.split("_vs_")[1]
        print n1
        print n2
        for k, v in dataframes.iteritems():
            if n1 in k:
                df1 = v  # selecting correct dataframe from the list of all dfs
                # print 'df 1 selected: ', len(df1)
            if n2 in k:
                df2 = v
    overlaps = pd.DataFrame(overlap)
    overlapingDF_row1 = overlaps['Sample1_row'].tolist()
    overlapingDF_row2 = overlaps['Sample2_row'].tolist()

    diff_df1 = pd.DataFrame()
    diff_df2 = pd.DataFrame()
    for index, value in df1.iterrows():
        #print index
        if index in overlapingDF_row1:
            continue
            #ind_overlap = overlapingDF_row1.index(index)
            #print overlaps['FC_value'][ind_overlap]
            #diff_df1 = diff_df1.append(value)
            #if overlaps['FC_value'][ind_overlap] >= 4 or overlaps['FC_value'][ind_overlap] <= 0.25:
                #diff_df1 = diff_df1.append(value)
        else:
            diff_df1 = diff_df1.append(value)
    for index, value in df2.iterrows():
        # print pos
        if index in overlapingDF_row2:
            continue
            #ind_overlap = overlapingDF_row2.index(index)
            #diff_df2 = diff_df2.append(value)
            #if index in set(overlapingDF_row2) and overlaps['FC_value'][ind_overlap] >= 4 or overlaps['FC_value'][ind_overlap] <= 0.25:
                #diff_df2 = diff_df2.append(value)
        else:
            diff_df2 = diff_df2.append(value)
    #unchanged_df = overlaps[overlaps['FC_value'] > 0.25]
    #unchanged_df1 = unchanged_df[unchanged_df['FC_value'] < 4]
    unchanged_df1 = overlaps
    diff_df1.to_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/'+ n1 + '_vs_' + n2 +'_Unique_'+ n1 + '.csv', sep=",", encoding='utf-8')
    diff_df2.to_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/Unique_'+ n1 + '_vs_' + n2 +'_Unique_'+ n2 + '.csv', sep=",", encoding='utf-8')
    #unchanged_df1.to_csv(
    #        '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/Unchanged_' + n1 + '_vs_' + n2 + '.csv', sep=",", encoding='utf-8')
    return diff_df1, diff_df2, unchanged_df1




















