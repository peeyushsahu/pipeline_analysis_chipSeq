__author__ = 'sahu'
import pandas as pd
import pysam
import timeit
from plotsAndseq import plots
import os

class PromoterEnrichment:
    '''
    To compute promoter specific enrichment for a ChIP-seq sample.
    '''
    def __init__(self, name, bampath, controlbam=None, gtf_transcript=None, path2save=''):
        self.name = name
        self.bampath = bampath
        self.controlbam = controlbam
        self.enrichmentdict = ''
        self.gtf_transcript = gtf_transcript
        self.path2save = path2save

    def heatmap_chip_enrichment(self):
        '''
        This will create a heatmap for all the transcript to model enrichment of tags around tss.
        :return:
        '''
        df = pd.read_csv('/ps/imt/f/Genomes/geneAnotations/gtf_transcript4annotation.db',
                                          header=0, sep='\t')
        print('Nos rows in df with chromosome variants:', len(df))
        df['chr'] = df['chr'].astype(str)
        df = df[df['chr'].str.len() < 4]
        print('Nos rows in df after removing chromosome variants:', len(df))
        df.index = range(len(df))
        self.gtf_transcript = df
        columns = ['chr', 'start', 'stop', 'strand', 'gene_name', 'transcript_id']
        heatmapdf = generate_heatmap(self.bampath, self.gtf_transcript, control_bampath=self.controlbam)
        big_df = plots.kmeans_clustering(heatmapdf, 9, 1000)
        print('Len of heatmap df', len(big_df))
        for col in columns[::-1]:
            #print(col, peak_df[col])
            big_df.insert(0, col, self.gtf_transcript[col])

        big_df.to_csv(os.path.join(self.path2save, self.name+'_heatmap_promoter_sub_R2_B6_RA_norm_50bp.tsv'), sep="\t", encoding='utf-8')


def generate_heatmap(bam_path, dataframe, control_bampath='', dist4middle=3000, steps=120, bpdist=50):
    '''
    Returns dataframe for tag count distribution for overlapping peaks within 3000bp (+,-) from summit.
    This function also considers the gene transcrition direction.
    This can also subtract background using control or IgG.
    :param bam_peak1:
    :param overlap_df:
    :return:
    '''
    startT = timeit.default_timer()
    sample_bam = pysam.Samfile(bam_path, "rb")
    if control_bampath is not '':
        print('Reading control bam file')
        control_bam = pysam.Samfile(control_bampath, "rb")
        totalreads_control = control_bam.mapped
        print('Total number of alignment in control:', totalreads_control)
    else:
        control_bam = None
    total_mapped = sample_bam.mapped
    print('Total number of alignment in sample:', total_mapped)
    peak_distribution_sample = pd.DataFrame()

    # check if genome of alignment (UCSC or ENSEMBL) bam
    try:
        sample_bam.count('9', 99181564, 99181974)
    except ValueError:
        print('Bam file is UCSC aligned converting coordinates accordingly...')
        dataframe['chr'] = 'chr'+dataframe['chr']
        pass

    overlap_df = dataframe[['chr', 'start', 'stop', 'strand', 'gene_name', 'transcript_id']]
    print('Process: Feature extraction from BAM started')

    for ind, row in overlap_df.iterrows():
        chr = str(row['chr'])
        orientation = row['strand']
        if orientation == '+':
            middle = row['start']
        if orientation == '-':
            middle = row['stop']

        start = middle - dist4middle  # Distance on one side of the peaks
        stop = start + bpdist
        list_sample1 = []
        # total_tags = int(bam_peak1.mapped) will get total no of mapped reads
        if start > 0:
            for i in range(0, steps):  # Please set based on distance on one side = s*distance/50
                seqcount = sample_bam.count(chr, start, stop)
                n_seqcount = (seqcount/total_mapped)*1000000
                # Calculating control tag count
                if control_bam is None:
                    list_sample1.append(seqcount)
                else:
                    controlcount = control_bam.count(chr, start, stop)
                    n_controlcount = (controlcount/totalreads_control)*1000000
                    list_sample1.append(seqcount-controlcount)
                start = stop
                stop = start + bpdist  # divide peaks into length of 50 bp
            if orientation == 1:  # Direction gene transcription
                # print 'Towards 5 prime'
                peak_distribution_sample = peak_distribution_sample.append(pd.Series(list_sample1),
                                                                           ignore_index=True)
            else:
                # print 'Towards 3 prime'
                peak_distribution_sample = peak_distribution_sample.append(pd.Series(list_sample1[::-1]),
                                                                           ignore_index=True)
        del middle, start, stop
    stop = timeit.default_timer()
    print('\nTime elapsed:' + str((stop - startT) / steps) + 'min')
    sample_bam.close()
    control_bam.close()
    return peak_distribution_sample

