__author__ = 'peeyush'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
import pandas as pd
import math
import numpy as np
import annotate.Annotate as Annotate
import overlap_analysis.differential_binding
from alignment import commons
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath


class PeaksAnalysis():
    """A dataframe for peak data analysis is created"""

    def __init__(self, peaks_df, con_name, dirPath):
        self.peaks = peaks_df
        self.name = con_name
        self.dirPath = dirPath
        self.GCount = None


    def genomic_regions(self):
        """
        This adds attribute to PeakAnalysis object for distribution of genomic regions
        :return:
        """
        GRcount = self.peaks['GenomicPosition TSS=1250 bp, upstream=5000 bp'].value_counts()
        GPcount = zip(GRcount.index, GRcount.values)
        self.GCount = GPcount
        with open(basepath + "/further_analysis/plots/G_distributions.txt", "a") as file:
            file.write(self.name+'\t'+str(GPcount)+'\n')
        #plotGenomicregions(GPcount, self.name)
        #stacked_plot_regions(GRcount.values, GPcount, self.name)

    def plot_factors(self, columnname):
        '''
        This function will plot a pie chart of occurrence for the unique elements in selected column.
        '''
        column = self.peaks[columnname].astype(basestring)
        factor = column.value_counts(sort=True)
        name = factor.index.tolist()
        y = factor.values.tolist()
        x = np.arange(len(y))
        #print '\n',x,'\n',y,'\n',name
        fig, ax = plt.subplots()
        plt.bar(x, y, color='coral')
        ax.set_xticks(x+0.5)
        ax.set_xticklabels(name, rotation=90)
        ax.set_xlim(0, len(name) + 0.5)
        plt.tight_layout()
        plt.title(self.name)
        plt.savefig(os.path.join(self.dirPath, columnname+'_'+self.name + '.png'))
        plt.clf()

def plotGenomicregions(GPcount, name):
    """
    :param GPcount: is a list of tuples [(region, size),....()]
    :return:
    """
    """ Now we produce some pie charts """
    gr = ['tss', 'intergenic', 'intron', 'exon', 'upstream']
    size = [0, 0, 0, 0, 0]
    for a, b in GPcount:
        if a == 'tss':
            size[0] = b
        if a == 'intergenic':
            size[1] = b
        if a == 'intron':
            size[2] = b
        if a == 'exon':
            size[3] = b
        if a == 'upstream':
            size[4] = b
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'cyan']
    explode = (0.1, 0, 0, 0, 0)  # only "explode" the 2nd slice
    plt.pie(size, explode=explode, labels=gr, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    #plt.legend(['tss', 'intergenic', 'intron', 'exon', 'upstream'], loc='upper left')
    plt.axis('equal')
    plt.savefig(basepath + '/further_analysis/plots/' + name + '.png')
    plt.clf()


def sumzip(*items):
    return [sum(values) for values in zip(*items)]


def stacke_plot_multiple(names_list, filtered_peaks, path):
    '''
    Creating Genomic Region stacked plot for multiple samples
    :param names_list:
    :param filtered_peaks:
    :return:
    '''
    tss = [7]
    intergenic = [46]
    intron = [42]
    exon = [2]
    upstream = [3]
    samples = ['Human']
    for name in names_list:
        samples.append(name.split(' ')[0])
        df = filtered_peaks.get(name)
        GRcount = df['GenomicPosition TSS=1250 bp, upstream=5000 bp'].value_counts()
        per = [100.0/sum(GRcount.values)*i for i in GRcount.values]
        GPcount = zip(GRcount.index, per)
        for G,C in GPcount:
            if G == 'tss':
                tss.append(C)
            if G == 'intergenic':
                intergenic.append(C)
            if G == 'intron':
                intron.append(C)
            if G == 'exon':
                exon.append(C)
            if G == 'upstream':
                upstream.append(C)

    N = len(names_list)+1
    ind = np.arange(N)    # the x locations for the groups
    width = 0.25       # the width of the bars: can also be len(x) sequence

    fig, ax = plt.subplots(figsize=(4, 6))
    p1 = plt.bar(ind, intergenic,   width, color='yellowgreen')
    p2 = plt.bar(ind, upstream, width, color='gold',
                 bottom=intergenic)
    p3 = plt.bar(ind, tss, width, color='lightskyblue',
                 bottom=sumzip(intergenic, upstream))
    p4 = plt.bar(ind, exon, width, color='lightcoral',
                 bottom=sumzip(intergenic, upstream, tss))
    p5 = plt.bar(ind, intron, width, color='cyan',
                 bottom=sumzip(intergenic, upstream, tss, exon))

    def autolabel(rects, gr_list):
        # attach some text labels
        i = 0
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., gr_list[i]+(height/2)-1, '%d'%int(height),
                    ha='center', va='bottom')
            gr_list[i] = gr_list[i] + height
            i += 1
        return gr_list

    gr_list = [0]*(len(names_list)+1)
    gr_list = autolabel(p1, gr_list)
    gr_list = autolabel(p2, gr_list)
    gr_list = autolabel(p3, gr_list)
    gr_list = autolabel(p4, gr_list)
    autolabel(p5, gr_list)

    plt.xlabel('Samples')
    plt.ylabel('% genomic distribution')
    plt.title('Genomic region ratio')
    plt.xticks(ind+width/2., samples)
    plt.ylim(0,100)
    plt.yticks(np.arange(0, 100, 100/5))
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('Intron', 'Exon', 'TSS', 'Upstream(-1500bp)', 'Intergenic(-1500 to 5000bp)') ,loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.savefig(os.path.join(path, 'multi_stacked' + ','.join(samples) + '.png'), bbox_inches='tight')
    #plt.clf()
    plt.close()


def peakTSSbinning(names, filtered_peaks, path):
    '''
    This plots a histogram for peak to TSS distance
    :param names:
    :param filtered_peaks:
    :param path:
    :return:
    '''
    import collections
    import matplotlib.pyplot as plt
    binSize = 100
    keys = range(-5000, 5000, binSize)
    #bins = collections.OrderedDict.fromkeys(keys)
    bins = dict(zip(keys, [0]*len(keys)))
    for k, v in filtered_peaks.iterrows():
        pos = v['Next Transcript tss distance']
        if (pos < 4900) and (pos > 0):
            key = int(math.ceil(pos / 100.0)) * 100
            bins[key] += 1
        if (pos > -5000) and (pos < 0):
            key = int(math.ceil(pos / 100.0)) * 100
            bins[key] += 1
        '''
        if pos < -5000:
            bins[-5000] += 1
        if pos > 4900:
            bins[4900] += 1
        '''
    #print(bins)
    keys.remove(0); keys.append(5000)
    bins = collections.OrderedDict(sorted(bins.items(), key=lambda t:t[0]))
    ## plotting barplot
    plt.subplots(figsize=(20, 6))
    plt.bar(range(len(bins)), bins.values(), align='center', color='#8A2BE2')
    plt.xticks(range(len(bins)), keys, rotation='vertical')
    plt.xlabel('TSS')
    plt.ylabel('Peak density')
    plt.title('Peak densities relative to TSS')
    plt.savefig(os.path.join(path, 'Peak_densities' + names + '.png'), bbox_inches='tight')
    #plt.clf()
    plt.close()

def OverlappingPeaks(dict_peaksdf, name, name1):

    """
    :param second_df: object of PeakAnalysis
    name: name+'vs'+name1
    :return: A dictionary of list of overlapping regions list(dict) and name
    """
    import timeit
    print 'Check point: Overlapping analysis'
    for k, v in dict_peaksdf.iteritems():
        if k == name:
            df1 = v
            #int name
        if k == name1:
            df2 = v
            #print 'size of df2: ', len(df2)
            #print name1
    print '\n',name,'vs',name1
    df1 = df1.peaks.sort(['chr'], ascending=True)
    df2 = df2.peaks.sort(['chr'], ascending=True)
    ### Method test PeakOverlaps
    start1 = timeit.default_timer()
    try:
        overlap_list = PeakOverlaps(df1, df2)
    except:
        print '\nWarning: Dataframe does not contain all the columns required for overlap, ' \
              'switching to minimal column requirement.'
        overlap_list = PeakOverlaps_concise(df1, df2)
    stop1 = timeit.default_timer()
    print "Time consumed by method PeakOverlaps:", stop1-start1, 'sec'
    overlap_dict = {name+'_vs_'+name1: overlap_list}
    ddf = pd.DataFrame(overlap_list)
    dirPath = os.path.join(basepath, 'further_analysis', 'overlap', name, '_vs_', name1)
    commons.ensure_path(dirPath)
    get_unique_peaks(df1, df2, name, name1, ddf, dirPath)
    ddf.to_csv(os.path.join(dirPath, name, '_vs_', name1, '.txt'), sep="\t", encoding='utf-8')
    return overlap_dict


def PeakOverlaps_concise(df1, df2):
    '''
    Used only when dataframe does not contain all the column needed.
    We will calculate overlap between peaks. An overlap is, if the two peaks share at-least one bp on genomic region.
    :param df1:
    :param df2:
    :return:
    '''

    df1['chr'] = df1['chr'].astype(str)
    df2['chr'] = df2['chr'].astype(str)
    df1_g = df1.groupby('chr')
    df2_g = df2.groupby('chr')
    num_overlap = 0
    overlap_list = []
    ind = 0
    j = 0
    for i in df1_g:
        #print '\nchr:', i[0]
        if i[0] in df2_g.groups:
            for count, row in df1_g.get_group(i[0]).iterrows():
                #print count
                j += 1
                if j == math.ceil(len(df1)/float(100)):
                    j = 0
                    ind += 1
                sys.stdout.write("\r%d%%" % ind)
                sys.stdout.flush()
                for count1, row1 in df2_g.get_group(i[0]).iterrows():
                    if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                          #
                          #
                          # it is an overlap 2
                          #
                          #
                        overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                        'Next transcript gene name':row['Next transcript gene name'], 'start1':row1['start'], 'stop1':row1['stop'], 'overlap':2, 'summit':row['summit'], 'summit1':row1['summit'],
                        'Repeat_name':row1['repeat'], 'repeat_family':row1['class/family']}
                        #overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                        #'Next transcript gene name':row['Next transcript gene name'], 'start1':row1['start'], 'stop1':row1['stop'],
                        #'Next transcript gene name1':row1['Next transcript gene name'], 'overlap':1, 'summit':row['summit'], 'summit1':row1['summit']}
                        overlap_list.append(overlaps)
                        num_overlap += 1
                        break
    return overlap_list

def PeakOverlaps(df1, df2):
    '''
    We will calculate overlap between peaks. An overlap is, if the two peaks share at-least one bp on genomic region.
    :param df1:
    :param df2:
    :return:
    '''

    df1['chr'] = df1['chr'].astype(str)
    df2['chr'] = df2['chr'].astype(str)
    df1_g = df1.groupby('chr')
    df2_g = df2.groupby('chr')
    num_overlap = 0
    overlap_list = []
    ind = 0
    j = 0
    for i in df1_g:
        #print '\nchr:',i[0]
        if i[0] in df2_g.groups:
            for count, row in df1_g.get_group(i[0]).iterrows():
                #print count
                j += 1
                if j == math.ceil(len(df1)/float(100)):
                    j = 0
                    ind += 1
                sys.stdout.write("\r%d%%" % ind)
                sys.stdout.flush()
                for count1, row1 in df2_g.get_group(i[0]).iterrows():
                    if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                          #
                          #
                          # it is an overlap 2
                          #
                          #
                        overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                        'Next transcript gene name':row['Next transcript gene name'], 'Next Transcript tss distance':row['Next Transcript tss distance'], 'start1':row1['start'], 'stop1':row1['stop'],
                        'Next transcript gene name1':row1['Next transcript gene name'], 'summit':row['summit'], 'summit1':row1['summit']}
                        overlap_list.append(overlaps)
                        num_overlap += 1
                        break
    return overlap_list


def get_unique_peaks(dataframe1, dataframe2, name, name1, overlapdf, dirpath):
    '''
    Write unique peaks for overlaping samples
    :param dataframe1:
    :param dataframe2:
    :param overlapdf:
    :param dirpath:
    :return:
    '''
    df1_overlap = list(overlapdf['Sample1_row'])
    df2_overlap = list(overlapdf['Sample2_row'])
    uni_df1 = dataframe1[~dataframe1.index.isin(df1_overlap)]
    uni_df2 = dataframe2[~dataframe2.index.isin(df2_overlap)]
    uni_df1.to_csv(os.path.join(dirpath, name, '_unique.txt'), sep='\t', index=None, header=True)
    uni_df2.to_csv(os.path.join(dirpath, name1, '_unique.txt'), sep='\t', index=None, header=True)


def non_overlapping_peaks(dataframe1, overlapDF):
    """
    Reterive non overlapping data from peak data
    :param overlapingDF:
    :param df1:
    :param df2:
    :return: tuple of non-overlapping data
    """
    print('Writing unique peaks....')
    print('Size of Df1:',len(dataframe1))
    print('Size of Df2:',len(overlapDF))
    uniqueDF = pd.DataFrame(columns=list(dataframe1.columns))
    dataframe1['chr'] = dataframe1['chr'].astype(str)
    overlapDF['chr'] = overlapDF['chr'].astype(str)
    df1_g = dataframe1.groupby('chr')
    df2_g = overlapDF.groupby('chr')
    for i in df1_g:
        #print '\nchr:', i[0]
        if i[0] in df2_g.groups:
            for k, v in df1_g.get_group(i[0]).iterrows():
                sys.stdout.write("\r%s%%" % i[0])
                sys.stdout.flush()
                in_overlap = False
                for ind, row in df2_g.get_group(i[0]).iterrows():
                    if (v['chr'] == row['chr']) & (v['start'] == row['start']):
                        #print('here')
                        in_overlap = True
                        break
                if not in_overlap:
                    uniqueDF = uniqueDF.append(v)
    uniqueDF = uniqueDF.rename(columns={'Next Gene name':'Next transcript gene name'})
    print('Size of Unique Df1:',len(uniqueDF))
    return uniqueDF


def significance_of_chipseq_overlap(overlap, peak_df1, peak_df2, iteration=10):
    from statistics import mean
    import sys
    '''
    A permutation method to access the significance of overlap.
    :return:
    '''
    n_peak_df1 = peak_df1[['chr', 'start', 'stop']]
    n_peak_df2 = peak_df2[['chr', 'start', 'stop']]
    overlaps = overlap
    n_peak_df1['chr'] = n_peak_df1['chr'].astype(str)
    n_peak_df2['chr'] = n_peak_df2['chr'].astype(str)
    df = n_peak_df1.append(n_peak_df2, ignore_index=True)
    datapoints = min(len(n_peak_df1), len(n_peak_df2))
    print n_peak_df1.shape, n_peak_df2.shape, df.shape, datapoints
    ## Number of iterations
    overlap_list = {}
    for ite in range(0, iteration):
        sys.stdout.write("\r%d%%" % ite)
        sys.stdout.flush()
        df1_g = df.sample(n=datapoints)
        df2_g = df.loc[~df.index.isin(df1_g.index)]
        #print len(df1_g), len(df2_g)
        df1_g = df1_g.groupby(df1_g['chr'])
        df2_g = df2_g.groupby(df2_g['chr'])
        num_overlap = 0
        for i in df1_g:
            #print len(i[1])
            if i[0] in df2_g.groups:
                #print len(df2_g.get_group(i[0]))
                for count, row in df1_g.get_group(i[0]).iterrows():
                    for count1, row1 in df2_g.get_group(i[0]).iterrows():
                        if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                            #print row['start'], row1['start']
                            num_overlap += 1
                            break
        #print 'Overlap for iteration', num_overlap
        overlap_list[ite] = num_overlap
    #calculate p-value for the iterations
    pval = len([elem for elem in overlap_list.values() if elem >= overlap])/iteration
    return overlap_list, pval