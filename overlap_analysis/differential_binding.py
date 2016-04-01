__author__ = 'peeyush'

import pandas as pd
import os
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath


class Overlaps():
    """
    This class object will hold filteredPeaks data and peaksList for calculation differential binding.
    """

    def __init__(self, overlappinglist, filter_peaks):
        self.samples_names = overlappinglist
        self.filter_peaks = filter_peaks


    def diffBinding(self, basepeakfile, outpath=None, genewide=False):
        '''
        This function will extract summit (+-500) peak data if peak length is >1000 from provided peaks.
        This dataframe can be used with DESeq for differential binding calculation.
        :return:
        '''
        import pysam
        import sys, math
        print "\nCheck point: diffBinding"
        sample_name = self.samples_names
        dataframes = self.filter_peaks
        if genewide:
            print("Gene-wide calculation is on....")
            longestTranscriptDB = '/ps/imt/f/Genomes/geneAnotations/longest_transcript_annotation.db'
            transcriptDB = pd.read_csv(longestTranscriptDB, header=0, sep='\t')
            df = pd.DataFrame()
            df = pd.concat([df, dataframes.get(basepeakfile)['Next transcript gene name']], axis=1)
        # print type(sample_name)
        #df = dataframes.get(sample_name[0]).iloc[:, 0:1]
        else:
            df = pd.DataFrame()
            df = pd.concat([df, dataframes.get(basepeakfile)['chr']], axis=1)
            df = pd.concat([df, dataframes.get(basepeakfile)['start']], axis=1)
            df = pd.concat([df, dataframes.get(basepeakfile)['stop']], axis=1)
            df = pd.concat([df, dataframes.get(basepeakfile)['GenomicPosition TSS=1250 bp, upstream=5000 bp']], axis=1)
            df = pd.concat([df, dataframes.get(basepeakfile)['Next transcript gene name']], axis=1)
            df = pd.concat([df, dataframes.get(basepeakfile)['Next transcript strand']], axis=1)
            df = pd.concat([df, dataframes.get(basepeakfile)['Next Transcript tss distance']], axis=1)
            df = pd.concat([df, dataframes.get(basepeakfile)['summit']], axis=1)
            df['new_start'] = 0
            df['new_stop'] = 0
        #print df.head()
        print df.dtypes
        for sample in sample_name:
            df[sample] = 0
            #print '\n'+sample
            sample_bam_path = getBam(sample) #.split(' vs ')[0]
            sample_bam = pysam.Samfile(sample_bam_path, "rb")
            total_reads = sample_bam.mapped
            for k, v in df.iterrows():
                ## Gene wide differential calculation
                if genewide:
                    gene_name = v['Next transcript gene name']
                    if gene_name in list(transcriptDB['gene_name']):    # checking coordinates of genes in the database
                        gene_ind = transcriptDB['gene_name'][transcriptDB['gene_name'] == gene_name].index[0]
                        #print transcriptDB.loc[gene_ind,:]
                        chr = transcriptDB.loc[gene_ind, 'chr']; start = transcriptDB.loc[gene_ind, 'start']
                        stop = transcriptDB.loc[gene_ind, 'stop']; strand = transcriptDB.loc[gene_ind, 'strand']
                        tags = sample_bam.count(chr, start, stop)
                        if tags == 0: tags = 1
                        #print tags, gene_name, start, stop
                        df.loc[k,'chr'] = chr
                        df.loc[k,'start'] = start
                        df.loc[k,'stop'] = stop
                        df.loc[k,'strand'] = strand
                        df.loc[k,sample] = tags
                        df.loc[k,sample+'_norm_millon'] = (float(tags)/total_reads)*10**6
                        df.loc[k,sample+'_length_norm'] = ((float(tags)/total_reads)*10**6)/((float(stop)-start)/100)

                else:
                    sys.stdout.write("\rNumber of peaks processed:%d" % k)
                    sys.stdout.flush()
                    #print v['start'], v['summit']
                    if v['stop']-v['start'] > 500:
                        chr = str(v['chr'])
                        summit = v['start']+v['summit']
                        tags = sample_bam.count(chr, summit-500, summit+500)
                        if tags == 0: tags = 1
                        df.loc[k,'new_start'] = summit-1000
                        df.loc[k,'new_stop'] = summit+1000
                        df.loc[k,sample] = tags
                        df.loc[k,sample+'_norm_millon'] = (float(tags)/total_reads)*10**6
                    else:
                        chr = str(v['chr'])
                        summit = v['start']+v['summit']
                        tags = sample_bam.count(chr, v['start'], v['stop'])
                        if tags == 0: tags = 1
                        df.loc[k,'new_start'] = v['start']
                        df.loc[k,'new_stop'] = v['stop']
                        df.loc[k,sample] = tags
                        df.loc[k,sample+'_norm_millon'] = (float(tags)/total_reads)*10**6
            sample_bam.close()
        ## This will calculate the logFC
        for ind, row in df.iterrows():
            control = sample_name[0]
            #print row
            for sample in sample_name[1:]:
                df.loc[ind, 'log2FC_'+control+'_norm_vs_'+sample+'_norm'] = math.log(row[sample+'_norm_millon']/row[control+'_norm_millon'],2)

        if outpath is not None:
            outpath = outpath
        else:
            outpath = basepath + '/further_analysis/differential/' + '_'.join(sample_name) + '.txt'
        df.to_csv(outpath, sep="\t", encoding='utf-8', ignore_index=True)

def getBam(name):
    '''
    This function takes the sample name and return the path of its associated bam file.
    :param name:
    :return:
    '''
    from os import listdir
    path = basepath +'/results/AlignedLane/'
    bam_list = listdir(basepath + '/results/AlignedLane')
    Dir = None
    file = None
    print name
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
        return os.path.join(Dir, file)




def group_DF(dataframe, factor):
    """
    This function will sort the dataframe and return reletive positions of factor (eg. chromosome, genomic region) in peak table.
    :param overlapsObj:
    :return:
    """
    if not factor in dataframe.columns:
        raise ValueError('Not valid column for df groupping.')
    print 'DF grouping based on: '+factor
    grouppedDF = dataframe.groupby(factor)
    for k, v in grouppedDF:
        print 'Cluster:', k, 'Size:', len(v)
    return grouppedDF

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


def modification4nearestgenes(dataframe, name, modification_list):
    '''
    This function will extract summit (+-500) peak data if peak length is >1000 from provided peaks.
    This dataframe can be used with DESeq for differential bound calculation.
    :return:
    '''
    import pysam
    import sys
    print "\nCheck point: diffBinding"
    sample_name = name
    dataframes = dataframe
    # print type(sample_name)
    #df = dataframes.get(sample_name[0]).iloc[:, 0:1]
    df = pd.DataFrame()

    ### create df from nearest genes
    print 'Reassembling dataframe'
    genewidpos = zip(dataframes['next5genes'], dataframes['position'])
    for i in genewidpos:
        gene = i[0].split(',')
        pos = i[1].split(',')
        for j in range(0,len(gene)):
            row = [gene[j], pos[j].split(':')[0], pos[j].split(':')[1], pos[j].split(':')[2], (int(pos[j].split(':')[2])-int(pos[j].split(':')[1]))]
            #print row
            df = df.append(pd.Series(row), ignore_index=True)
    df.columns = ['gene', 'chr', 'start', 'stop', 'length']

    ### Extract tag count for the specific region of peak
    for sample in modification_list:
        print '\n'+sample+' sample being processed.'
        df[sample] = 0
        sample_bam_path = getBam(sample) #.split(' vs ')[0]
        sample_bam = pysam.Samfile(sample_bam_path, "rb")
        if 'pol' in sample or 'K4me3' in sample:
            for k, v in df.iterrows():
                sys.stdout.write("\rNumber of peaks processed:%d" % k)
                sys.stdout.flush()
                chr = str(v['chr'])
                tags = sample_bam.count(chr, int(v['start'])-500, int(v['start'])+500)
                df.loc[k, sample] = float(tags)/1000
        else:
            for k, v in df.iterrows():
                #print v['chr'], v['start'], v['stop']
                sys.stdout.write("\rNumber of peaks processed:%d" % k)
                sys.stdout.flush()
                chr = str(v['chr'])
                tags = sample_bam.count(chr, int(v['start']), int(v['stop']))
                norm_tags = float(tags)/(int(v['stop']) - int(v['start']))
                df.loc[k, sample] = norm_tags
        sample_bam.close()
    df.to_csv(
            basepath + '/further_analysis/differential/' + sample_name + '_'.join(modification_list) + '.csv',
            sep=",", encoding='utf-8', ignore_index=True)




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
            basepath + '/further_analysis/differential/'+ n1 + '_vs_' + n2 +'_Unique_'+ n1 + '.csv', sep=",", encoding='utf-8')
    diff_df2.to_csv(
            basepath + '/further_analysis/differential/Unique_'+ n1 + '_vs_' + n2 +'_Unique_'+ n2 + '.csv', sep=",", encoding='utf-8')
    #unchanged_df1.to_csv(
    #        '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/Unchanged_' + n1 + '_vs_' + n2 + '.csv', sep=",", encoding='utf-8')
    return diff_df1, diff_df2, unchanged_df1




















