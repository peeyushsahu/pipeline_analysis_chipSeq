import gc
from overlap_analysis.differential_binding import getBam, factor_seperate


__author__ = 'peeyush'


def make_dir(bam_order, region='All'):
    import os
    # print 'Directory_for_result: ' + '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/'+folder
    path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlapping_plots/' + bam_order + '/' + region + '/'
    if not os.path.exists(path):
        os.makedirs(path)
    return path


#@profile
def GR_heatmaps_DF_for_peaks(bam_name_list, peak_df, region=None, sort=False, sort_column=None):
    '''
    Suggestion: Please do not use more then 3 samples at a time.
    This function will take a list of bam files path and create a heatmap of genomic region for the provided peak dataset.
    :param bam_name_list: A list cantaining names of bam files to be vizualized eg. bam_name_list = ['PRMT6_2_seq6', 'H3K4me3_seq2', 'Sample_K9me3']
    :param peak_df: A datframe containing peak information
    :param region: Which region to plot (eg. tss, intron, exon, intergenic or None)
    :param sort: True if dataframe should be sorted
    :param sort_column: If sorted=True; give the column name to be sorted
    :return: A dataframe; Column: additive length of genomic regions for all the bam files, Rows: Number of peaks defined in the peak dataframe
    '''
    import pandas as pd
    import gc
    big_df = pd.DataFrame()

    region = region.strip()
    sort = sort
    peak_df = peak_df
    sort_column = sort_column
    if region != 'all':
        peak_df = peak_df[peak_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == region]

    if region == 'tss':  # Reduce peaks based on their distance from TSS
        print 'Selecting peaks only within +-300bp'
        peak_df = peak_df[peak_df['Next Transcript tss distance'] < 300]
        peak_df = peak_df[peak_df['Next Transcript tss distance'] > -300]
        if len(peak_df) == 0:
            raise ValueError('selected region does not contain any peaks')
    print region + ' found in dataframe: ', len(peak_df)
    #print peak_df.head()
    peak_df.index = range(0, len(peak_df))
    if sort:
        colnames = peak_df.columns.tolist()
        indices1 = [i for i, s in enumerate(colnames) if sort_column in s]
        for i in indices1:
            if "RA" not in colnames[i] and "RA" not in sort_column and "norm" not in colnames[i]:
                condition = colnames[i]
                print 'Sorted on column: ' + condition
                peak_df = peak_df.sort(condition, ascending=True)
                sort_column = condition
                break
            elif "RA" in colnames[i] and "RA" in sort_column and "norm" not in colnames[i]:
                condition = colnames[i]
                print 'Sorted on column: ' + condition
                peak_df = peak_df.sort(condition, ascending=True)
                sort_column = condition
                break

    #print peak_df.head()
    for v in bam_name_list:
        df = overlapping_peaks_distribution(v, peak_df)
        df = scale_dataframe(df)   # scaling of dataframe
        big_df = pd.concat([big_df, df], axis=1)
    #big_df = overlapping_peaks_distribution_improved(bam_name_list, peak_df)
    big_df.columns = range(0, big_df.shape[1])
    bam_order = ','.join(bam_name_list)
    path = make_dir(bam_order, region + str(len(peak_df)))
    #big_df.to_csv(path+bam_order+region+'.csv', sep=",", encoding='utf-8', ignore_index=True)
    #divide_peaks_in_strength(big_df, bam_order, path)         # plotting line plots
    try:
        big_df = kmeans_clustering(big_df, 9, 100)  # performing k-means clustering
        dict_of_df = factor_seperate(big_df, 'cluster')  # divide df in smaller dfs basis in clustering
        line_plot_peak_distribution(dict_of_df, bam_order, path)  # plotting individual clusters
        if len(bam_name_list) == 4:
            print 'No. of sample to plot: 4'
            plot_clustered_peaks_4_four_samples(dict_of_df, bam_order,
                                                path)  # plotting cluster for different bam in overlapping plot
        if len(bam_name_list) == 3:
            plot_clustered_peaks_4_three_samples(dict_of_df, bam_order, path)
        if len(bam_name_list) == 2:
            plot_clustered_peaks_4_two_samples(dict_of_df, bam_order, path)
    except:
        print 'Dataframe can not be clustered, scipy error.'

    ### adding columns to heatmap df
    try:
        if sort:big_df.insert(0, sort_column, peak_df[sort_column])
        big_df.insert(0, 'GenomicPosition TSS=1250 bp, upstream=5000 bp', peak_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'])
        big_df.insert(0, 'Next transcript gene name', peak_df['Next transcript gene name'])
        big_df.insert(0, 'Next transcript strand', peak_df['Next transcript strand'])
        big_df.insert(0, 'summit', peak_df['summit'])
        big_df.insert(0, 'stop', peak_df['stop'])
        big_df.insert(0, 'start', peak_df['start'])
        big_df.insert(0, 'chr', peak_df['chr'])
    except:
        big_df.insert(0, 'Next transcript strand', peak_df['Next transcript strand'])
        big_df.insert(0, 'summit', peak_df['summit'])
        big_df.insert(0, 'stop', peak_df['stop'])
        big_df.insert(0, 'start', peak_df['start'])
        big_df.insert(0, 'chr', peak_df['chr'])
    #print big_df.head()
    big_df.to_csv(path + bam_order + region + '.csv', sep=",", encoding='utf-8') #, ignore_index=True
    gc.collect()

def kmeans_clustering(df, nClus, iter):
    '''
    This will perform clustering on a dataframe
    :param df: DataFrame
    :param nClus: No. of clusters
    :param iter: Max no. of iterations
    :return: Dataframe attached with cluster information column
    '''
    import scipy.cluster.vq as cluster
    import pandas as pd

    print '\nProcess: Clustering of DataFrame'
    df_fin = df
    df_mat = df.as_matrix()
    res = cluster.kmeans2(df_mat, nClus, iter)
    df_fin.insert(0, 'cluster', pd.Series(res[1]))
    return df_fin

#@profile
def overlapping_peaks_distribution(bam_name, overlap_df):
    '''
    Returns dataframe for tag count distribution for overlapping peaks within 500bp (+,-) from summit.
    This function also considers the gene transcrition direction.
    :param bam_peak1:
    :param overlap_df:
    :return:
    '''
    import pysam
    import pandas as pd
    import sys
    import timeit
    startT = timeit.default_timer()
    bam_path = getBam(bam_name)
    sample_bam = pysam.Samfile(bam_path, "rb")
    peak_distribution_sample = pd.DataFrame()
    overlap_df = overlap_df[['chr', 'start', 'stop', 'Next transcript strand', 'summit']]
    print 'Process: Feature extraction from BAM started'
    count = 1
    summit_found = ''
    for ind in overlap_df.T.to_dict().values():
        sys.stdout.write("\rFeature extraction for peak:%d" % count)
        sys.stdout.flush()
        count += 1
        chr = str(ind['chr'])
        if 'chr' in chr:
            chr = chr[3:]
        orientation = ind['Next transcript strand']
        middle = ind['start'] + ind['summit']
        start = middle - 3000   # Distance on one side of the peaks
        stop = start + 50
        list_sample = []
        # total_tags = int(bam_peak1.mapped) will get total no of mapped reads
        if start > 0:
            for i in range(0, 120):  # Please set based on distance on one side = s*distance/50
                seqcount = sample_bam.count(chr, start, stop)
                list_sample.append(seqcount)
                start = stop
                stop = start + 50  # divide peaks into length of 50 bp
            if orientation == 1:  # Direction gene transcription
                # print 'Towards 5 prime'
                peak_distribution_sample = peak_distribution_sample.append(pd.Series(list_sample), ignore_index=True)
            else:
                # print 'Towards 3 prime'
                peak_distribution_sample = peak_distribution_sample.append(pd.Series(list_sample[::-1]), ignore_index=True)
        del list_sample, middle, start, stop
    stop = timeit.default_timer()
    print summit_found
    print '\nTime elapsed:'+str((stop-startT)/60)+'min'
    sample_bam.close()
    return peak_distribution_sample


def divide_peaks_in_strength(df, name, path):
    '''
    This will divide peaks on the basis of strength (eg. low, mid-mid, mid-high, high) and plots it using def plot_divide_peaks_in_strength()
    :param df:
    :param name:
    :return:
    '''
    lenght = len(df) / 4
    start = 0
    end = lenght
    dict_list = {}
    for i in range(0, 4):
        sub_df = df[start:end]
        dict_list[i] = sub_df
        start = end - 1
        end += lenght
    # print dict_list
    plot_divide_peaks_in_strength(dict_list, name, path)
    return dict_list


def plot_divide_peaks_in_strength(dict_list, name, path):
    '''
    Plots result of divide_peaks_in_strength.
    :param dict_list:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    # print 'shape', dict_list[0].shape
    x = np.array(range(dict_list[0].shape[1] / 2 * -50, dict_list[0].shape[1] / 2 * 50, 50))
    # print x
    s = np.array(dict_list[0].sum(axis=0))
    l = np.array(dict_list[1].sum(axis=0))
    m = np.array(dict_list[2].sum(axis=0))
    h = np.array(dict_list[3].sum(axis=0))

    plt.ylim(0, max(max(h), max(m), max(s), max(l)) + 100)
    plt.xlabel('Binding profile')
    plt.ylabel('Normalized tag density')
    plt.title('Genomic distribution of peaks')
    plt.gca().set_color_cycle(['mediumorchid', 'coral', 'r', 'dodgerblue'])  #'mediumorchid', 'coral',

    xnew = np.linspace(x.min(), x.max(), 300)
    smooth = spline(x, s, xnew)
    plt.plot(xnew, smooth, linewidth=3)

    xnew1 = np.linspace(x.min(), x.max(), 300)
    smooth1 = spline(x, l, xnew1)
    plt.plot(xnew1, smooth1, linewidth=3)

    xnew2 = np.linspace(x.min(), x.max(), 300)
    smooth2 = spline(x, m, xnew2)
    plt.plot(xnew2, smooth2, linewidth=3)

    xnew3 = np.linspace(x.min(), x.max(), 300)
    smooth3 = spline(x, h, xnew3)
    plt.plot(xnew3, smooth3, linewidth=3)

    plt.legend(['Low', 'Medium-Low', 'Medium-high', 'High'], loc='upper left')  #'Low', 'Medium',
    #plt.show()
    plt.savefig(path + name + '.png')
    plt.clf()


def line_plot_peak_distribution(dict_of_df, name, path):
    '''
    This will plot clusters individually
    :param dict_of_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np

    for Cluster, df in dict_of_df.iteritems():
        plt.figure(figsize=(10.5, 8.5))
        df = df.drop('cluster', axis=1)
        size = df.shape
        # x_axis = []
        # for i in range(0, size[1]/80):
        #    x_axis.extend(range(-2000,+2000,50))
        #x = np.array(x_axis)

        x = np.array(range(df.shape[1] / 2 * -50, df.shape[1] / 2 * 50, 50))
        #print len(x)
        s = np.array(df.sum(axis=0))
        #print len(s)
        plt.ylim(0, max(s) + 50)
        plt.xlabel('Distribution of cluster: ' + str(Cluster))
        plt.ylabel('Binding profile cluster ' + str(Cluster))
        plt.title('Cluster associated peak distribution with datapoints:' + str(size[0]))
        plt.gca().set_color_cycle(['r'])
        xnew = np.linspace(x.min(), x.max(), 300)
        smooth = spline(x, s, xnew)
        plt.plot(xnew, smooth, linewidth=3)  #marker='o'
        #plt.legend([names[1], names[3]], loc='upper left')
        #plt.show()
        plt.savefig(path + name + '_cluster:' + str(Cluster) + '.png')
        plt.savefig(path + name + '_cluster:' + str(Cluster) + '.svg')
        plt.clf()


def plot_clustered_peaks_4_two_samples(dict_df, name, path):
    '''
    Plots result of divided_cluster_peaks_in_strength.
    :param dict_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    import scipy as sp
    # print 'shape', dict_list[0].shape
    for k, v in dict_df.iteritems():
        size = v.shape
        size_one_df = size[1]/2
        df1 = v.iloc[:, 1:size_one_df+1]
        df2 = v.iloc[:, size_one_df+1:size[1]]
        #x = np.array(range(-3000, 3000, 50))
        x = np.array(range(df1.shape[1] / 2 * -50, df1.shape[1] / 2 * 50, 50))
        # print x
        s = np.array(df1.sum(axis=0))/float(len(df1))
        #s = np.subtract(s, sp.median(s))           ### Normalizing with median
        #s = np.subtract(s, np.percentile(s, 25))    ### Normalizing with 25 percentile
        s = np.subtract(s, min(s))

        l = np.array(df2.sum(axis=0))/float(len(df2))
        #l = np.subtract(l, sp.median(l))
        #l = np.subtract(l, np.percentile(l, 25))
        l = np.subtract(l, min(l))

        plt.ylim(0, max(max(s), max(l)) + 10)
        plt.xlabel('Binding profile cluster' + str(k))
        plt.ylabel('Combined tag density')
        plt.title('Genomic distribution of peaks with datapoints: ' + str(size[0]))
        plt.gca().set_color_cycle(['mediumorchid', 'r', 'dodgerblue'])  #'mediumorchid', 'coral',

        xnew = np.linspace(x.min(), x.max(), 300)
        smooth = spline(x, s, xnew)
        plt.plot(xnew, smooth, linewidth=3)

        xnew1 = np.linspace(x.min(), x.max(), 300)
        smooth1 = spline(x, l, xnew1)
        plt.plot(xnew1, smooth1, linewidth=3)

        sname = name.split(',')
        plt.legend([sname[0], sname[1]], loc='upper left')  #'Low', 'Medium',
        #plt.show()
        plt.savefig(path + 'overlap_' + name + '_cluster:' + str(k) + '.png')
        plt.savefig(path + 'overlap_' + name + '_cluster:' + str(k) + '.svg')
        plt.clf()


def plot_clustered_peaks_4_three_samples(dict_df, name, path):
    '''
    Plots result of divided_cluster_peaks_in_strength.
    :param dict_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    import scipy as sp
    # print 'shape', dict_list[0].shape
    for k, v in dict_df.iteritems():
        size = v.shape
        df1 = v.iloc[:, 1:121]
        df2 = v.iloc[:, 121:241]
        df3 = v.iloc[:, 241:361]
        x = np.array(range(-3000, 3000, 50))
        # print x
        s = np.array(df1.sum(axis=0))/float(len(df1))
        #s = np.subtract(s, sp.median(s))           ### Normalizing with median
        #s = np.subtract(s, np.percentile(s, 25))    ### Normalizing with 25 percentile
        s = np.subtract(s, min(s))

        l = np.array(df2.sum(axis=0))/float(len(df2))
        #l = np.subtract(l, sp.median(l))
        #l = np.subtract(l, np.percentile(l, 25))
        l = np.subtract(l, min(l))

        m = np.array(df3.sum(axis=0))/float(len(df3))
        #m = np.subtract(m, sp.median(m))
        #m = np.subtract(m, np.percentile(m, 25))
        m = np.subtract(m, min(m))

        plt.ylim(0, max(max(m), max(s), max(l)) + 10)
        plt.xlabel('Binding profile cluster' + str(k))
        plt.ylabel('combined tag density / no. of peaks')
        plt.title('Genomic distribution of peaks with datapoints: ' + str(size[0]))
        plt.gca().set_color_cycle(['mediumorchid', 'r', 'dodgerblue'])  #'mediumorchid', 'coral',

        xnew = np.linspace(x.min(), x.max(), 300)
        smooth = spline(x, s, xnew)
        plt.plot(xnew, smooth, linewidth=3)

        xnew1 = np.linspace(x.min(), x.max(), 300)
        smooth1 = spline(x, l, xnew1)
        plt.plot(xnew1, smooth1, linewidth=3)

        xnew2 = np.linspace(x.min(), x.max(), 300)
        smooth2 = spline(x, m, xnew2)
        plt.plot(xnew2, smooth2, linewidth=3)

        sname = name.split(',')
        plt.legend([sname[0], sname[1], sname[2]], loc='upper left')  #'Low', 'Medium',
        #plt.show()
        plt.savefig(path + 'overlap_' + name + '_cluster:' + str(k) + '.png')
        plt.savefig(path + 'overlap_' + name + '_cluster:' + str(k) + '.svg')
        plt.clf()


def plot_clustered_peaks_4_four_samples(dict_df, name, path):
    '''
    Plots result of divided_cluster_peaks_in_strength.
    :param dict_df:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import spline
    import numpy as np
    import scipy as sp
    # print 'shape', dict_list[0].shape
    for k, v in dict_df.iteritems():
        size = v.shape
        df1 = v.iloc[:, 1:121]
        df2 = v.iloc[:, 121:241]
        df3 = v.iloc[:, 241:361]
        df4 = v.iloc[:, 361:481]
        x = np.array(range(-3000, 3000, 50))
        # print x
        s = np.array(df1.sum(axis=0))/float(len(df1))
        #s = np.subtract(s, sp.median(s))           ### Normalizing with median
        #s = np.subtract(s, np.percentile(s, 25))    ### Normalizing with 25 percentile
        s = np.subtract(s, min(s))

        l = np.array(df2.sum(axis=0))/float(len(df2))
        #l = np.subtract(l, sp.median(l))
        #l = np.subtract(l, np.percentile(l, 25))
        l = np.subtract(l, min(l))

        m = np.array(df3.sum(axis=0))/float(len(df3))
        #m = np.subtract(m, sp.median(m))
        #m = np.subtract(m, np.percentile(m, 25))
        m = np.subtract(m, min(m))

        n = np.array(df4.sum(axis=0))/float(len(df4))
        n = np.subtract(n, min(n))

        plt.ylim(0, max(max(m), max(s), max(l), max(n)) + 10)
        plt.xlabel('Binding profile cluster' + str(k))
        plt.ylabel('Combined tag density')
        plt.title('Genomic distribution of peaks with datapoints: ' + str(size[0]))
        plt.gca().set_color_cycle(['mediumorchid', 'r', 'dodgerblue', 'green'])  #'mediumorchid', 'coral',

        xnew = np.linspace(x.min(), x.max(), 300)
        smooth = spline(x, s, xnew)
        plt.plot(xnew, smooth, linewidth=3)

        xnew1 = np.linspace(x.min(), x.max(), 300)
        smooth1 = spline(x, l, xnew1)
        plt.plot(xnew1, smooth1, linewidth=3)

        xnew2 = np.linspace(x.min(), x.max(), 300)
        smooth2 = spline(x, m, xnew2)
        plt.plot(xnew2, smooth2, linewidth=3)

        xnew3 = np.linspace(x.min(), x.max(), 300)
        smooth3 = spline(x, n, xnew3)
        plt.plot(xnew3, smooth3, linewidth=3)

        sname = name.split(',')
        plt.legend([sname[0], sname[1], sname[2], sname[3]], loc='upper left')  #'Low', 'Medium',
        #plt.show()
        plt.savefig(path + 'overlap_' + name + '_cluster:' + str(k) + '.png')
        plt.savefig(path + 'overlap_' + name + '_cluster:' + str(k) + '.svg')
        plt.clf()


def scale(val, src, dst):
    '''
    This returns scaled value
    :param val: value to be scaled
    :param src: max and min of values to be scaled
    :param dst: range to scale
    :return:
    '''
    return (((val - src[0]) * (dst[1] - dst[0])) / (src[1] - src[0])) + dst[0]


def scale_list(list, Min, Max, oMin, oMax):
    '''
    This takes list and returns lst of scaled values
    :param list:
    :param Min:
    :param Max:
    :return:
    '''
    old_list = list
    # old_max = max(old_list)
    # old_min = min(old_list)
    old_max = oMax
    old_min = oMin
    new_max = Max
    new_min = Min
    new_list = []
    for i in old_list:
        new_list.append(scale(i, (old_min, old_max), (new_min, new_max)))
    return new_list


def scale_dataframe(df):
    '''
    This will scale a Pandas.DataFrame
    :param df: DataFrame
    :return: Scaled DataFrame
    '''
    import pandas as pd

    old_max = int(max(df.max(axis=1, numeric_only=True)))
    old_min = 0
    new_max = 100
    new_min = 0
    list_of_rows = []
    print 'Process: scaling of dataframe'
    print 'Max value in df:',old_max
    if old_max > 50:  # Scale values only when the highest in dataframe is > 50
        for r, v in df.iterrows():
            rows = []
            for val in v:
                # print val
                if not isinstance(val, str):
                    # print val
                    rows.append(scale(val, (old_min, old_max), (new_min, new_max)))
            list_of_rows.append(rows)
        scaled_df = pd.DataFrame(data=list_of_rows)
    else: scaled_df = df
    return scaled_df


def plot_overlapping_peaks(peak_distribution_sample1, peak_distribution_sample2, name):
    '''
    This function plots the results of "overlapping_peaks_distribution()" into a line plot.
    :param peak_distribution_sample1:
    :param peak_distribution_sample2:
    :param name:
    :return:
    '''
    import matplotlib.pyplot as plt

    list_a = []
    list_b = []
    for i in range(0, peak_distribution_sample1.shape[1]):
        list_a.append(int(sum(peak_distribution_sample1[i])))
        list_b.append(int(sum(peak_distribution_sample2[i])))

    scale_a = scale_list(list_a, 0, 100)
    scale_b = scale_list(list_b, 0, 100)
    # print len(scale_a), len(scale_b)
    x = range(-1200, 1200, 50)
    # print x
    plt.ylim(0, 120)
    plt.xlabel('TSS')
    plt.ylabel('Normalized peak height')
    plt.title(name)
    plt.gca().set_color_cycle(['dodgerblue', 'mediumorchid'])
    plt.plot(x, scale_a, linewidth=3)
    plt.plot(x, scale_b, linewidth=3)
    plt.legend(['PRMT6', 'H3K4'], loc='upper left')
    # plt.show()
    plt.savefig('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/' + name + '.png')
    plt.clf()


def overlapping_peak_intensity(peakdf1, name1, peakdf2, name2, overlappingdf, overlap_name):
    # ## Reference bins created for original peak data
    import collections as cl

    colist1 = peakdf1.columns.values.tolist()
    colist2 = peakdf2.columns.values.tolist()

    indices1 = [i for i, s in enumerate(colist1) if name1 in s]
    indices2 = [i for i, s in enumerate(colist2) if name2 in s]

    df1_tags = sorted(peakdf1.ix[:, indices1[1]], reverse=False)
    df2_tags = sorted(peakdf2.ix[:, indices2[1]], reverse=False)

    df1_bin = []
    df1_length = len(peakdf1) / 10
    start = 0
    end = df1_length
    for i in range(0, 10):
        df1_bin.append(sum(df1_tags[start:end]) / df1_length)
        print 'df1', sum(df1_tags[start:end]) / df1_length
        start = end
        end = start + df1_length

    df2_bin = []
    df2_length = len(peakdf2) / 10
    start = 0
    end = df2_length
    for i in range(0, 10):
        df2_bin.append(sum(df2_tags[start:end]) / df2_length)
        print 'df2', sum(df2_tags[start:end]) / df2_length
        start = end
        end = start + df2_length

    # ## Overlapping bins created for overlapping peaks
    import re

    overlap_name_list = re.split(" vs |_vs_", overlap_name)
    o_colist = overlappingdf.columns.values.tolist()
    print overlap_name_list
    o_indices1 = [s for i, s in enumerate(o_colist) if 'norm_' + overlap_name_list[0] in s]
    o_indices2 = [s for i, s in enumerate(o_colist) if 'norm_' + overlap_name_list[2] in s]
    print o_indices1

    overlappingdf = overlappingdf.sort(o_indices1[0], ascending=True)
    o_df1_tags = overlappingdf[o_indices1]
    o_df2_tags = overlappingdf[o_indices2]

    ol_df1_bin = []
    df_length = len(overlappingdf) / 10
    start = 0
    end = df_length
    for i in range(0, 10):
        ol_df1_bin.append(o_df1_tags[start:end].sum()[0] / df_length)
        print 'df1', o_df1_tags[start:end].sum()[0] / df_length
        start = end
        end = start + df_length

    ol_df2_bin = []
    start = 0
    end = df_length
    for i in range(0, 10):
        ol_df2_bin.append(o_df2_tags[start:end].sum()[0] / df_length)
        print 'df2', o_df2_tags[start:end].sum()[0] / df_length
        start = end
        end = start + df_length

    List = cl.OrderedDict()
    List[name1] = df1_bin
    List[name2] = df2_bin
    List['Overlap_' + overlap_name_list[0]] = ol_df1_bin
    List['Overlap_' + overlap_name_list[2]] = ol_df2_bin

    # line_plot_overlapping_peak_intensity(List)
    line_plot_overlapping_peak_intensity_2(List)
    return List


def line_plot_overlapping_peak_intensity(dict_of_bins):
    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    import matplotlib.pyplot as plt

    if 1:
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.75)

        par1 = host.twinx()
        par2 = host.twinx()
        # par3 = host.twinx()

        offset = 40
        new_fixed_axis = par2.get_grid_helper().new_fixed_axis
        par2.axis["right"] = new_fixed_axis(loc="right",
                                            axes=par2,
                                            offset=(offset, 0))
        # new_fixed_axis = par3.get_grid_helper().new_fixed_axis
        #par3.axis["right"] = new_fixed_axis(loc="right",
        #                                    axes=par3,
        #                                    offset=(2 * offset, 0))

        par2.axis["right"].toggle(all=True)
        #par3.axis["right"].toggle(all=True)

        List = dict_of_bins.values()
        names = dict_of_bins.keys()
        x_range = range(0, len(List[0]))

        host.set_xlim(0, len(List[0]))
        host.set_ylim(0, int(max(List[1])) + 10)

        host.set_xlabel("Clustered peaks")
        host.set_ylabel(names[1])
        par1.set_ylabel(names[2])
        par2.set_ylabel(names[3])
        #par3.set_ylabel(names[3])

        p1, = host.plot(x_range, List[1], label=names[1], marker='o')
        p2, = par1.plot(x_range, List[2], label=names[2], marker='o')
        p3, = par2.plot(x_range, List[3], label=names[3], marker='o')
        #p4, = par3.plot(x_range, List[3], label=names[3], marker='o')

        par1.set_ylim(0, int(max(List[2])) + 10)
        par2.set_ylim(0, int(max(List[3])) + 10)
        #par3.set_ylim(0, int(max(List[3])) + 10)

        host.legend(loc='upper left')

        host.axis["left"].label.set_color(p1.get_color())
        par1.axis["right"].label.set_color(p2.get_color())
        par2.axis["right"].label.set_color(p3.get_color())
        #par3.axis["right"].label.set_color(p4.get_color())

        plt.draw()
        # plt.show()
        plt.savefig(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/overlapping_peak_intensity_' +
            names[0] + names[1] + '.png')
        plt.clf()


def line_plot_overlapping_peak_intensity_2(dict_of_bins):
    import matplotlib.pyplot as plt

    dict_list = dict_of_bins.values()
    names = dict_of_bins.keys()
    # s = dict_list[0]
    # x = range(len(dict_list[1]))
    l = dict_list[1]
    m = dict_list[2]
    h = dict_list[3]
    plt.ylim(0, max(l) + 10)
    plt.xlabel('Normalized PRMT6 peaks')
    plt.ylabel('Corresponding H3k4me3 peaks')
    plt.title('H3K4 wrt PRMT6 ascending')
    plt.gca().set_color_cycle(['dodgerblue', 'r'])  #'mediumorchid', 'coral',

    #plt.plot(m, l, linewidth=3, marker='o')

    #xnew1 = np.linspace(x.min(),x.max(),300)
    #smooth1 = spline(x,l,xnew1)
    #plt.plot(xnew1, smooth1, linewidth=3)

    #xnew2 = np.linspace(x.min(),x.max(),300)
    #smooth2 = spline(x,m,xnew2)
    #plt.plot(xnew2, smooth2, linewidth=3)

    plt.plot(m, h, linewidth=3, marker='o')

    plt.legend([names[1], names[3]], loc='upper left')  #'Low', 'Medium',
    #plt.show()
    plt.savefig('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/' + names[2] + '_' + names[
        3] + '.png')
    plt.clf()


def peak_position_dataframe(peak_df, name):
    '''
    This function will take dataframe and extract the peak position onto a matrix from -2000 to 2000 on genome
    :param peak_df:
    :return:
    '''
    import pandas as pd
    import math

    geneNames = []
    positionDf = pd.DataFrame()
    for k, v in peak_df.iterrows():
        # print v['Next Transcript tss distance']
        if v['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'tss':
            positionList = [0] * 40
            position = int(math.ceil(v['Next Transcript tss distance'] / 100))
            # print position
            posOnList = 20 + position
            if posOnList >= 0 and posOnList <= 40:
                #print posOnList
                geneNames.append(v['Next transcript gene name'])
                positionList[posOnList - 1] = 1
                positionDf = positionDf.append(pd.Series(positionList), ignore_index=True)
        # positionDf['Next transcript gene name'] = pd.Series(geneNames)
        positionDf = positionDf.set_index(pd.Series(geneNames))
        #plot_heatmap_4_peaks_position(positionDf)
    print 'Frequency of peak positions\n', positionDf.sum()
    positionDf.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/' + name + '.csv',
                      sep=",", encoding='utf-8', ignore_index=True)
    return positionDf

