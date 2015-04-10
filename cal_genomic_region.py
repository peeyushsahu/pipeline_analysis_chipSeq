__author__ = 'peeyush'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import math

class PeaksAnalysis():
    """A dataframe for peak data analysis is created"""

    def __init__(self, peaks_df, con_name):
        self.peaks = peaks_df
        self.name = con_name
        self.GCount = None


def genomic_regions(self):

    """
    This adds attribute to PeakAnalysis object for distribution of genomic regions
    :return:
    """

    GPcount = self.peaks['GenomicPosition TSS=1250 bp, upstream=5000 bp'].value_counts()
    GPcount = zip(GPcount.index, GPcount.values)
    self.GCount = GPcount
    with open("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/plots/G_distributions.txt", "a") as file:
        file.write(self.name+'\t'+str(GPcount)+'\n')
    plotGenomicregions(GPcount, self.name)

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
    plt.savefig('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/plots/' + name + '.svg')
    plt.clf()


def PeakOverlaps(df1, df2):
    '''
    We will calculate overlap between peaks. An overlap is, if the two peaks share at-least one bp on genomic region.
    :param df1:
    :param df2:
    :return:
    '''
    num_overlap = 0
    overlap_list = []
    i = 0
    j = 0
    for count, row in df1.iterrows():
        #print count
        j += 1
        if j == math.ceil(len(df1)/float(100)):
            j = 0
            i += 1
        sys.stdout.write("\r%d%%" % i)
        sys.stdout.flush()
        for count1, row1 in df2.iterrows():
            #print type(str(row['chr'])) == type(str(row1['chr']))
            if str(row['chr']) == str(row1['chr']):
                #print row['start'], row1['start'], row['stop'], row1['stop']
                if (row['start'] >= row1['start']) and (row['stop'] <= row1['stop']):
                    #
                    #
                    # it is an complete overlap 1
                    #
                    #
                    overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                'Next transcript gene name':row['Next transcript gene name'], 'Next Transcript tss distance':row['Next Transcript tss distance'], 'start1':row1['start'], 'stop1':row1['stop'],
                'Next transcript gene name1':row1['Next transcript gene name'], 'overlap':1, 'summit':row['summit'], 'summit1':row1['summit']}
                    overlap_list.append(overlaps)
                    num_overlap += 1
                    #print num_overlap
                    break


                if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                  #
                  #
                  # it is an overlap 2
                  #
                  #
                    overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                'Next transcript gene name':row['Next transcript gene name'], 'Next Transcript tss distance':row['Next Transcript tss distance'], 'start1':row1['start'], 'stop1':row1['stop'],
                'Next transcript gene name1':row1['Next transcript gene name'], 'overlap':2, 'summit':row['summit'], 'summit1':row1['summit']}
                    overlap_list.append(overlaps)
                    num_overlap += 1
                    break

                if max(row['start'], row1['start'])-1000 < min(row['stop'], row1['stop']):
                  #
                  #
                  # it is an overlap 3
                  #
                  #
                    overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                'Next transcript gene name':row['Next transcript gene name'], 'Next Transcript tss distance':row['Next Transcript tss distance'], 'start1':row1['start'], 'stop1':row1['stop'],
                'Next transcript gene name1':row1['Next transcript gene name'], 'overlap':3, 'summit':row['summit'], 'summit1':row1['summit']}
                    overlap_list.append(overlaps)
                    num_overlap += 1
                    break

    return overlap_list

def OverlappingPeaks(self, name, name1):

    """
    :param second_df: object of PeakAnalysis
    name: name+'vs'+name1
    :return: A dictionary of list of overlapping regions list(dict) and name
    """
    print 'Check point: Overlapping analysis.'
    for k, v in self.iteritems():
        if k == name:
            df1 = v
            #int name
        if k == name1:
            df2 = v
            #print 'size of df2: ', len(df2)
            #print name1
    #print type(df1)
    df1 = df1.peaks.sort(['chr'], ascending=True)
    #print 'size of df1: ', len(df1)
    df2 = df2.peaks.sort(['chr'], ascending=True)
    #print 'size of df2: ', len(df2)
    #print df1.columns.values.tolist()
    overlap_list = PeakOverlaps(df1, df2)
    overlap_dict = {name+'_vs_'+name1: overlap_list}
    ddf = pd.DataFrame(overlap_list)
    #ddf.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/check_'+name+'_vs_'+name1+'.csv', sep=",", encoding='utf-8')
    return overlap_dict


def overlappingPeaksLess(overlapping_samples, filtered_peak_data, name1, name2):
    """
    This function is same as the overlappingPeaks; only this can not be used for differential binding analysis.
    :return:
    """
    print 'Check point: overlappingPeaksLess'
    print type(overlapping_samples)
    df1 = []
    df2 = []
    for k, v in overlapping_samples.iteritems():
        #print 'This is key: '+k
        if k == name1:
            df1 = pd.DataFrame(v)
            #print name1
            #print type(df1)
        if k == name2:
            df2 = pd.DataFrame(v)
            #print name2

    if len(df1) < 1:
        #print 'df1 from peak_data'
        df1 = filtered_peak_data.get(name1)
    if len(df2) < 1:
        #print 'df2 from peak data'
        df2 = filtered_peak_data.get(name2)
    print type(df1), type(df2)
    overlap_list = PeakOverlaps(df1, df2)
    overlap_dict = {name1+'_vs_'+name2: overlap_list}
    ddf = pd.DataFrame(overlap_list)
    ddf.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/'+name1+'_vs_'+name2+'.csv', sep=",", encoding='utf-8')
    return overlap_dict