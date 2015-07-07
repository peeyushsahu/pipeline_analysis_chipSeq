__author__ = 'peeyush'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import math
import numpy as np

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
        plt.tight_layout()
        plt.title(self.name)
        plt.savefig('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/plots/' + columnname +'_'+self.name + '.png')
        plt.clf()

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
    ddf.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/'+name+'_vs_'+name1+'.csv', sep=",", encoding='utf-8')
    return overlap_dict

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
    plt.savefig('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/plots/' + name + '.png')
    plt.clf()


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
        print '\nchr:', i[0]
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
                    if (row['start'] >= row1['start']) and (row['stop'] <= row1['stop']):
                          #
                          #
                          # it is an complete overlap 1
                          #
                          #
                        overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                        'Next transcript gene name':row['Next transcript gene name'], 'start1':row1['start'], 'stop1':row1['stop'], 'overlap':1, 'summit':row['summit'], 'summit1':row1['summit'],
                        'Repeat_name':row1['repeat'], 'repeat_family':row1['class/family']}
                        overlap_list.append(overlaps)
                        num_overlap += 1
                        #print overlaps
                        break

                    if max(row['start'], row1['start']) < min(row['stop'], row1['stop']):
                          #
                          #
                          # it is an overlap 2
                          #
                          #
                        overlaps = {'Next transcript strand':row['Next transcript strand'],'Sample1_row':count, 'Sample2_row':count1, 'chr':row['chr'], 'start':row['start'], 'stop':row['stop'], 'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                        'Next transcript gene name':row['Next transcript gene name'], 'start1':row1['start'], 'stop1':row1['stop'], 'overlap':2, 'summit':row['summit'], 'summit1':row1['summit'],
                        'Repeat_name':row1['repeat'], 'repeat_family':row1['class/family']}
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
        print '\nchr:',i[0]
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
                        #print overlaps
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
    return overlap_list