__author__ = 'peeyush'

import pandas as pd
import rpy2


class Overlaps():
    """
    This class object will hold overlapping data for differential binding analysis
    """

    def __init__(self, overlappinglist, peak_data, filter_peaks):
        self.overlaps = overlappinglist
        self.peak_data = peak_data
        self.filter_peaks = filter_peaks
        self.chr_pos = None
        self.db_overlap = None  #holds object from diffBinding calculations


def diffBinding(self):
    print "\nCheck point: diffBinding"
    iter_overlaps = self.overlaps
    dataframes = self.peak_data
    # print type(iter_overlaps)
    for name, overlap in iter_overlaps.iteritems():
        print name.split("_vs_")
        n1 = name.split("_vs_")[0]  # split to get the sample name
        n2 = name.split("_vs_")[1]
        for k, v in dataframes.iteritems():
            if n1 in k:
                df1 = v  # selecting correct dataframe from the list of all dfs
                # print 'df 1 selected: ', len(df1)
            if n2 in k:
                df2 = v
                # print 'df 2 selected', len(df2)
        # print type(df1)
        # print n1.split(' ')
        nc1 = n1.split(' ')[0]  # further name splitting for selection of rowname in dataframe
        nc2 = n2.split(' ')[0]
        # print nc1
        colist1 = df1.columns.values.tolist()
        colist2 = df2.columns.values.tolist()
        #print colist1
        indices1 = [i for i, s in enumerate(colist1) if nc1 in s]
        indices2 = [i for i, s in enumerate(colist2) if nc2 in s]
        #print enumerate(colist1)
        #print indices1
        for row in overlap:
            #print row
            row1 = row.get('Sample1_row')
            #print row1
            row2 = row.get('Sample2_row')
            #print row2
            norm_val1 = df1.iloc[row1][indices1[1]]
            norm_val2 = df2.iloc[row2][indices2[1]]
            FC = norm_val1 / norm_val2  #sample one / sample two
            row['FC_value'] = FC
            row['Tag_count_norm_'+nc1] = norm_val1
            row['Tag_count_norm_'+nc2] = norm_val2
            row['Raw_Tag_count_'+nc1] = df1.iloc[row1][indices1[0]]
            row['Raw_Tag_count_'+nc2] = df2.iloc[row2][indices2[0]]
            #print FC
        self.db_overlaps = {name: overlap}
        overlap1 = pd.DataFrame(overlap)
        overlap1.to_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/' + n1 + '_vs_' + n2 + '.csv',
            sep=",", encoding='utf-8', ignore_index=True)

def diffBinding_p6(self):
    '''
    This will join the overlapping peak, e.g. min(start, start1) and max(stop, stop1).
    And will consider the bam files for tag-counts in the region.
    :return:
    '''
    import pysam
    print "\nCheck point: diffBinding"
    iter_overlaps = self.overlaps
    #dataframes = self.peak_data

    for name, overlap in iter_overlaps.iteritems():
        print name.split("_vs_")
        n1 = name.split("_vs_")[0]  # split to get the sample name
        n2 = name.split("_vs_")[1]
        nc1 = n1.split(' ')[0]  # further name splitting for selection of rowname in dataframe
        nc2 = n2.split(' ')[0]

        sample1_bam = pysam.Samfile(getBam(nc1), 'rb')
        sample2_bam = pysam.Samfile(getBam(nc2), 'rb')
        total_tags_s1 = int(sample1_bam.mapped)
        total_tags_s2 = int(sample2_bam.mapped)     # get sequencing depth
        print nc1+': ', total_tags_s1, ', '+nc2+': ', total_tags_s2
        for row in overlap:
            #print row2
            chr = str(row['chr'])
            #print chr
            start = int(min(row['start'], row['start1']))
            #print start
            stop = int(max(row['stop'], row['stop1']))
            #print stop
            tags = sample1_bam.count(chr, start, stop)
            tags1 = sample2_bam.count(chr, start, stop)
            #print tags, tags1
            ntags = (float(tags)/total_tags_s1)*10000000
            ntags1 = (float(tags1)/total_tags_s2)*10000000
            FC = float(tags) / tags1  #sample one / sample two
            #print FC
            row['FC_value'] = FC
            row['Tag_count_norm_'+nc1] = tags
            row['Tag_count_norm_'+nc2] = tags1
            row['Raw_Tag_count_'+nc1] = ntags
            row['Raw_Tag_count_'+nc2] = ntags1
            row['new_start'] = start
            row['new_stop'] = stop
            row['lenght'] = stop - start
            #print FC
        self.db_overlaps = {name: overlap}
        overlap1 = pd.DataFrame(overlap)
        sample1_bam.close()
        sample2_bam.close()
        overlap1.to_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/' + n1 + '_vs_' + n2 + '.csv',
            )

def getBam(name):
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
                        print 'Bam file selected: '+j
            if 'RA' not in name and 'RA' not in i:
                Dir = path+i
                #print Dir
                for j in listdir(Dir):
                    if j.endswith('.bam'):
                        file = j
                        print 'Bam file selected: '+j
    if file is None:
        raise KeyError('Bam file cannot be found for '+name)
    else:
        return Dir+'/'+file




def factor_seperate(dataframe, factor):
    """
    This function will return reletive positions of chromosomes in peak table
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


    # def getDf4OverlappingPeaks(overlappingDF, df1):
    """
    Retrieves individual peaks DF for overlapping peaks for PRMT6 and H3K4 peak analysis
    :return: DataFrame
    """
    # df1 = df1.set_index(df1['Unnamed: 0'])
    #overlapingDF_row1 = overlappingDF['Sample1_row']
    #df = df1[overlappingDF]



















