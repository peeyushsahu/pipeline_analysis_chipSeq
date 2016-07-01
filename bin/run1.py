import os
from pandas import read_csv
from alignment import retroT_analysis
from annotate.Annotate import next6genes_annotator
from overlap_analysis import cal_genomic_region, differential_binding, filterPeaks
from overlap_analysis.differential_binding import modification4nearestgenes
from plotsAndseq import seqOperations
import plotsAndseq.plots as plots
from plotsAndseq.plots import GR_heatmaps_DF_for_peaks
from plotsAndseq.seqOperations import density_based_motif_comparision
import pandas as pd
#import pdb; pdb.set_trace()
import gc
from stat_tests import permutation

__author__ = 'peeyush'
#time.strftime("%d/%m/%y")
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath

folders = ["overlap",
           "differential",
           "filtered",
           "plots",
           "seq4motif",
           "valley_peaks",
           "CpG",
           "density_based_motif",
           "combined"]
for folder in folders:
    #print 'Directory_for_result: ' + '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/'+folder
    path = basepath + '/further_analysis/'+folder
    if not os.path.exists(path):
        os.makedirs(path)
print 'Output folder created'

sample_name = [
    #'H3R2ame2_E9 vs IgG_E.9 filtered',
    #'H3R2me2a_E9_RA vs IgG_E9_RA filtered'
]

'''
['H3K4me3_B5.1 vs IgG_B5.1 filtered',
 'Sample_K27ac vs IgG_seq6 filtered',
 'PRMT6_E9_commer vs IgG_E.9 filtered',
 'H3K27me3_E9_RA vs IgG_E9_RA filtered',
 'Sample_K27ac_10_wt vs IgG_seq6_RA filtered',
 'H3R2me2a_B6.2_RA vs IgG_B6_RA filtered',
 'PRMT6_seq5 vs IgG_seq4 filtered',
 'PRMT6_seq6_RA vs IgG_seq6_RA filtered',
 'H3K27ac_E9_RA vs IgG_E9_RA filtered',
 'H3K36me3_B6.2_RA vs IgG_B6_RA filtered',
 'H3K4me3_B6.2_RA vs IgG_B6_RA filtered',
 'H3R2ame2_E9 vs IgG_E.9 MACS2 filtered',
 'H3K4me3_seq2 vs IgG_seq2 filtered',
 'Sample_EZH1 vs IgG_seq6 filtered',
 'H3K27me3_B6 vs IgG_B6.2 filtered',
 'PRMT6_B6_3ul vs IgG_B6.2 filtered',
 'H3K4me1_B6 vs IgG_B6.2 filtered',
 'H3K4me3_seq2 vs Sample_PIS filtered',
 'Sample_EZH2_RA vs IgG_seq6_RA filtered',
 'Sample_pol2_RA vs IgG_seq6_RA filtered',
 'PRMT6_KO_B5.1 vs IgG_B5.1 filtered',
 'Sample_K27ac_RA vs IgG_seq6_RA filtered',
 'Sample_K4me1 vs Sample_PIS filtered',
 'H3K27me3_E9 vs IgG_E.9 filtered',
 'Sample_K27me1_RA vs IgG_seq6_RA filtered',
 'H3K36me3_E9_RA vs IgG_E9_RA filtered',
 'Sample_K9me3_RA vs IgG_seq6_RA filtered',
 'H3K36me3_B6.2 vs IgG_B6.2 filtered',
 'H3K4me3_E9 vs IgG_E.9 filtered',
 'H3K27ac_B6.2 vs IgG_B6.2 filtered',
 'Sample_K4me1_RA vs Sample_PIS_RA filtered',
 'Sample_18F3 vs Sample_8C9 filtered',
 'H3K4me1_E9_RA vs IgG_E9_RA filtered',
 'PRMT6_1_seq1 vs IgG_seq2 filtered',
 'H3R2me2a_B6.2 vs IgG_B6.2 filtered',
 'PRMT6_KO_10.8 vs IgG_10.8 filtered',
 'Sample_K4me3_10_wt vs IgG_seq6_RA filtered',
 'Sample_PRMT6_2_10_wt vs IgG_seq6_RA filtered',
 'H3K27me3_B6_RA vs IgG_B6_RA filtered',
 'YY1_seq3 vs IgG_seq4 filtered',
 'Sample_EZH2 vs IgG_seq6 filtered',
 'PRMT6_B6_commer vs IgG_B6.2 filtered',
 'Sample_K36me3_RA vs IgG_seq6_RA filtered',
 'H3K4me1_B6_RA vs IgG_B6_RA filtered',
 'H3K27ac_B6_RA vs IgG_B6_RA filtered',
 'H3R2me2a_E9_RA vs IgG_E9_RA filtered',
 'H3K27ac_B5.1 vs IgG_B5.1 filtered',
 'H3K4me1_E9 vs IgG_E.9 filtered',
 'H3K4me3_seq2_RA vs Sample_PIS_RA filtered',
 'H3K27ac_E9 vs IgG_E.9 filtered',
 'H3K4me3_E9_RA vs IgG_E9_RA filtered',
 'Sample_PRMT6_3_0_wt vs IgG_seq6 filtered',
 'H3K36me3_E9 vs IgG_E.9 filtered',
 'PRMT6_seq5_RA vs IgG_seq4_RA filtered',
 'Sample_K27me3 vs IgG_seq6 filtered',
 'YY1_seq3_RA vs IgG_seq4_RA filtered',
 'PRMT6_seq4_RA vs IgG_seq4_RA filtered',
 'Sample_EZH1_RA vs IgG_seq6_RA filtered',
 'Sample_K9me3 vs IgG_seq6 filtered',
 'Sample_K27me3_RA vs IgG_seq6_RA filtered',
 'H3R2ame2_E9 vs IgG_E.9 filtered',
 'PRMT6_E9_3ul vs IgG_E.9 filtered',
 'PRMT6_seq6 vs IgG_seq6 filtered',
 'Sample_pol2 vs IgG_seq6 filtered',
 'Sample_K36me3 vs IgG_seq6 filtered',
 'PRMT6_seq4 vs IgG_seq4 filtered',
 'Sample_18F3_RA vs IgG_seq6_RA filtered',
 'PRMT6_KO_B6.2 vs IgG_B6.2 filtered',
 'Sample_K27me1 vs IgG_seq6 filtered',
 'H3R2ame2_B5.1 vs IgG_B5.1 filtered',
 'PRMT6_E9_3ul vs PRMT6_B6_3ul filtered',
 'Sample_PRMT6_3_RA vs IgG_seq6_RA filtered',
 'PRMT6_KO_E.9 vs PRMT6_KO_B6.2 filtered',
 'PRMT6_KO_E.9 vs IgG_E.9 filtered']

'''


# Here import peak called data in a list....
#filterPeaks.rawpeaks_in_allsamples()
file = open(basepath + '/further_analysis/filtered/filteredPeaksCount.txt', 'w')
file.write('Sample lane\t Background lane\tNo. of raw peaks\tNo. of filtered peaks\n')
Filter = True
peakAnalysis_df = {}
for name in sample_name:
    df = read_csv(basepath + '/csv/' + name + '.csv', header=0, sep='\t')
    df = df.rename(columns={'Next Gene name':'Next transcript gene name'})
    filtered_peak_data, dirPath = filterPeaks.filterpeaks(df, name, filtering=Filter)
    GR_analysis = cal_genomic_region.PeaksAnalysis(peaks_df=filtered_peak_data, con_name=name, dirPath=dirPath)
    GR_analysis.genomic_regions()
    peakAnalysis_df[name] = GR_analysis
    GR_analysis.plot_factors('Next Gene biotype')
    GR_analysis.plot_factors('chr')
    cal_genomic_region.stacke_plot_multiple([name], {name:filtered_peak_data}, dirPath)
    cal_genomic_region.peakTSSbinning(name, filtered_peak_data, dirPath)
print "Number of sample are being analysed: ", peakAnalysis_df.__len__()
print "Filtering peaks."
file.close()

'''
df = read_csv(
            basepath + '/further_analysis/PRMT6_KO_analysis/improved_PRMT6_E9_B6_B5_all_diff.txt',
            header=0, sep='\t')
df = df[df['log2FC_PRMT6_KO_E.9_norm_vs_PRMT6_KO_B6.2_norm'] < -0.8]
filtered_peak_data['PRMT6_peaks_improved'] = df
## Plot stacked plot for selected samples

cal_genomic_region.stacke_plot_multiple(['PRMT6_peaks_improved']
                                        , filtered_peak_data)

#cal_genomic_region.stacke_plot_multiple(['H3K4me3_seq2 vs IgG_seq2 filtered', 'H3K4me3_RA_seq2 vs IgG_RA_seq2 filtered']
#                                        , filtered_peak_data)
#cal_genomic_region.stacke_plot_multiple(['Sample_18F3 vs Sample_8C9 filtered', 'Sample_18F3_RA vs IgG_RA_seq6 filtered']
#                                        , filtered_peak_data)
'''

### Performs differential binding calulation from full sample
'''
sample = ['H3R2ame2_E9', 'H3R2me2a_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2']
peak_df_name = 'H3R2ame2_E9 vs IgG_E.9 filtered'
filtered_peak = {peak_df_name: peakAnalysis_df[peak_df_name].peaks}
diffbind = differential_binding.Overlaps(sample, filtered_peak)
diffbind.diffBinding(peak_df_name, outpath=basepath + '/further_analysis/PRMT6_KO_analysis/H3R2me2/H3R2me2a_2700/H3R2me2a_2700_E9_B6_wid_K4_diff.txt')
'''
### Diff.binding for nearest genes
'''
sample = ['Sample_K36me3', 'Sample_K36me3_RA',
          'Sample_pol-2', 'Sample_pol-2_RA',
          'H3K4me3_seq2', 'H3K4me3_RA_seq2',
          'Sample_K27me3', 'Sample_K27me3_RA']
peak_df = read_csv(basepath + '/further_analysis/differential/diffPeaks_P6_nearest_5_genes.csv',
    header=0, sep=',')

modification4nearestgenes(peak_df, 'prmt6_nearest5genes', sample)
'''

### Diff. Binding calculation from altered sample (external)
'''
# Perform overlap and perform diff binding
#overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, 'H3K4me3_E9 vs IgG_E.9 filtered', 'H3R2ame2_E9 vs IgG_E.9 filtered')

#'H3R2ame2_E9', 'H3K36me3_E9', 'H3K36me3_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2', 'H3K27ac_E9','H3K27ac_B6.2','PRMT6_E9_commer', 'PRMT6_B6_commer','PRMT6_KO_E.9', 'PRMT6_KO_B6.2','H3K4me1_E9','H3K4me1_B6','H3K27me3_E9','H3K27me3_B6'
path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/Compare -ATRA&+ATRA'
multiple_df = ['H3R2me2+-RA_K4me3_+RA_expr']

for df in multiple_df:
    peak_df = read_csv(os.path.join(path, df+'.txt'), header=0, sep='\t')
    peak_df = peak_df[(peak_df['log2FC'] <= -1) & (peak_df['pvalue'] <= 0.05)]
    print(peak_df.shape)
    #peak_df = peak_df[:10]
    sample = ['H3K27me3_E9', 'H3K27me3_B6', 'H3K27me3_E9_RA', 'H3K27me3_B6_RA']

    # random sampling of df
    #peak_df = differential_binding.random_sampleing_df(peak_df, 1459)

    peak_df = peak_df.rename(columns={'Next Gene name':'Next transcript gene name'})
    print peak_df.shape
    #peak_df = peak_df.drop_duplicates(['Next transcript gene name'])
    #print peak_df.shape
    filtered_peak = {'loaded_sample': peak_df}
    diffbind = differential_binding.Overlaps(sample, filtered_peak)
    outpath = basepath + '/further_analysis/H3R2me2a_analysis'
    diffbind.diffBinding('loaded_sample', outpath=outpath+'/H3R2me2_+-RA_H3K4me3_down_expr+RA_K27me3.txt', genewide=True, use_second_start=False, from_sample=2, twoK=False, highest=False) #, genewide=True
'''

### Transposone analysis
'''
trans = retroT_analysis.Transposones()
trans.read_repeat_pos('/ps/imt/e/RepBase20.05.fasta/rtro_human/transposone_db.csv')
trans.get_transposition_position(['SVA_D'])
region = ['all'] #'all', 'tss', 'exon', 'intron', 'intergenic', 'upstream'
bam_list = [['PRMT6_2_seq6', 'H3K4me3_seq2', 'Sample_K27ac', 'Sample_K4me1', 'Sample_K9me3']]
for bams in bam_list:
    for i in region:
        GR_heatmaps_DF_for_peaks(bams, trans.repeat_pos_db.get('SVA_D'), region=i,
                                 sort=False, sort_column='H3K4me3_RA_seq2')
'''
### Calculate overlap of Diff peaks with transposon
'''
diff_peaks = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlapping_plots/PRMT6_2_seq6,Sample_pol-2,Sample_K36me3,H3K4me3_seq2,Sample_K27me1,Sample_K27ac,Sample_K4me1,Sample_K9me3,Sample_K27me3,Sample_18F3/all27744'
                      '/PRMT6_2_seq6,Sample_pol-2,Sample_K36me3,H3K4me3_seq2,Sample_K27me1,Sample_K27ac,Sample_K4me1,Sample_K9me3,Sample_K27me3,Sample_18F3all.csv',
    header=0, sep=',')
diff_peaks = diff_peaks[(diff_peaks['cluster'] == 3) | (diff_peaks['cluster'] == 4) | (diff_peaks['cluster'] == 6) | (diff_peaks['cluster'] == 7)]

TE_df = read_csv('/ps/imt/e/RepBase20.05.fasta/rtro_human/transposone_db.csv',
    header=0, sep=',')
TE_df = TE_df[TE_df['class/family'] == 'Simple_repeat' | TE_df['class/family'] == 'Satellite/telo' | TE_df['class/family'] == 'Low_complexity']

#peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlapping_plots/Sample_K27ac,Sample_K4me1,H3K4me3_seq2/all22096/clustered_df_0,2,3.csv',
#                   header=0, sep=',')

### If DF is from R change column names ('' ','.')
if '.' in diff_peaks.columns[6]:
    from string import maketrans
    cols = diff_peaks.columns
    new_cols = []
    for col in cols:
        new_cols.append(col.translate(maketrans('.',' ')))
    diff_peaks.columns = new_cols


TE_df['chr'] = map(lambda x: x[3:], TE_df['chr'])
filtered_peak_data['TE_db'] = TE_df
filtered_peak_data['Total_PRMT6_TSS_3,4,6,7'] = diff_peaks

peakAnalysis_df = {}
for k, v in filtered_peak_data.iteritems():
    name = k
    df = v
    GR_analysis = cal_genomic_region.PeaksAnalysis(df, name)
    peakAnalysis_df[name] = GR_analysis

sample_name2 = ['TE_db']
sample_name1 = ['Total_PRMT6_TSS_3,4,6,7']

if len(sample_name1) != len(sample_name2):
    raise ValueError("Unequal sample list for comparison")
else:
    for i in range(0, len(sample_name1)):
        overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, sample_name1[i], sample_name2[i])
'''

### Compares multiple ChIP-Seq profile using peaks (heatmap) from on sample
'''
region = ['all'] #'all', 'tss', 'exon', 'intron', 'intergenic', 'upstream'
bam_list = [['H3R2ame2_E9', 'H3R2me2a_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2']]
for bams in bam_list:
    for i in region:
        df = peakAnalysis_df['H3R2ame2_E9 vs IgG_E.9 filtered'].peaks
        GR_heatmaps_DF_for_peaks(bams, df, region=i,
                                 sort=False, sort_column='H3R2ame2_E9', scale_df=False, normalized=True, strength_divide=False)
'''

### Combine peaks of two df
'''
df1 = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/H3R2me2a_E9_RA vs IgG_E9_RA filtered.csv', sep='\t', header=0, index_col=None)
df2 = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/H3R2me2a_B6.2_RA vs IgG_B6_RA filtered.csv', sep='\t', header=0, index_col=None)
peaks_df_n = cal_genomic_region.get_combined_peaks(df1,df2)
#sampleName = 'combined_H3R2me2a_E9+B6'
peaks_df_n.to_csv(basepath+'/further_analysis/combined/combined_H3R2me2a_E9+B6_RA.txt', sep='\t', header=True)
'''


### Gene-wide chip profile for broad histone marks
'''
bam_list = ['H3K36me3_E9', 'H3K36me3_B6.2']
sample_name = 'H3R2ame2_E9 vs IgG_E.9 filtered_vs_H3K4me3_E9 vs IgG_E.9 filtered'
peak_df = read_csv(basepath + '/further_analysis/H3R2me2a_analysis/H3R2ame2_E9 vs IgG_E.9 filtered_vs_H3K4me3_E9 vs IgG_E.9 filtered/'+sample_name+'.txt', sep='\t', header=0)
#peak_df = peak_df[:200]
plots.grHeatmap4wholeGene(peak_df, bam_list, sample_name)
'''

### Comapre ChIP-Seq profile from altered sample (external)

# Perform overlap
#overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, 'H3R2ame2_E9 vs IgG_E.9 filtered', 'H3K4me3_E9 vs IgG_E.9 filtered')
listGroups = ['H3R2ame2_E9 vs IgG_E.9 filtered']#,'chip_H3K4me1_up','chip_H3K4me3_down','chip_H3K4me3_up','chip_H3K27ac_down','chip_H3K27ac_up']

for sampleName in listGroups:
    peaks_df_n = read_csv(basepath + '/further_analysis/filtered/H3R2ame2_E9 vs IgG_E.9 filtered/'+sampleName+'.txt', sep='\t', header=0)
    #peaks_df_n = overlapping_res #overlapping_res
    #peaks_df_n = peaks_df_n[peaks_df_n['log2FC_PRMT6_E9_commer_norm_vs_PRMT6_B6_commer_norm'] <= -0.8]
    print 'Dim of DF', peaks_df_n.shape
    bam_list = [['H3R2ame2_E9', 'H3R2me2a_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2', 'H3K27me3_E9', 'H3K27me3_B6']]
    #'H3K4me3_E9', 'H3K4me3_B6.2', 'H3K4me3_E9_RA', 'H3K4me3_B6.2_RA','H3K27ac_E9', 'H3K27ac_B6.2', 'H3K27ac_E9_RA', 'H3K27ac_B6_RA', 'H3K4me1_E9_RA', 'H3K4me1_B6_RA'
    #'Sample_18F3_RA', 'Sample_18F3', 'H3R2ame2_E9','H3R2me2a_B6.2','H3K27ac_E9','H3K27ac_B6.2','H3K4me1_E9','H3K4me1_B6','H3K27me3_E9','H3K27me3_B6',
    #'H3K27me3_E9', 'H3K27me3_B6', 'H3R2me2a_E9_RA', 'H3R2me2a_B6.2_RA', 'H3R2ame2_E9', 'H3R2me2a_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2','H3K27ac_E9','H3K27ac_B6.2','H3K4me1_E9', 'H3K4me1_B6', 'PRMT6_KO_E.9', 'PRMT6_KO_B6.2', 'PRMT6_KO_B5.1'

    # If DF is from R change column names ('' ','.')
    peaks_df_n = peaks_df_n.rename(columns={'Next Gene name':'Next transcript gene name'})
    #cal_genomic_region.peakTSSbinning(sampleName, peaks_df_n, path=basepath + '/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/PRMT6_new+old_peaks')

    #peaks_df_n = peaks_df_n[peaks_df_n['cluster'].isin([4,8,12,16,17,19,20])]
    if '.' in peaks_df_n.columns[6]:
        from string import maketrans
        cols = peaks_df_n.columns
        new_cols = []
        for col in cols:
            new_cols.append(col.translate(maketrans('.',' ')))
        peaks_df_n.columns = new_cols

    for List in bam_list:
        region = ['all'] #'all', 'tss', 'exon', 'intron', 'intergenic', 'upstream'
        for i in region:
            GR_heatmaps_DF_for_peaks(List, peaks_df_n, region=i, sort=False, sort_column='PRMT6_KO_E.9_norm_millon', scale_df=False,
                                     sample_name=sampleName, strength_divide=False)
gc.collect()


### permutation test for difference in median/mean
'''
name = ['H3K4me3_E9','H3K4me3_B6.2','H3K27ac_E9','H3K27ac_B6.2']
comparisions = {-0.64:['H3K27ac_E9', 'H3K27ac_B6.2'],
                0.05:['H3K4me3_E9', 'H3K4me3_B6.2']}
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3K4me3_all_+K27ac+-RA.txt', sep='\t', header=0)
print(peak_df.head())
outdir = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis'
permutation.permutation_test4peakdensity(peak_df, name, comparisions, n=310, niter=10000, outdir=outdir, sname='H3K4me3_down')
'''



### Density based motif analysis
'''
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/overlapping/PRMT6_KO_E.9,PRMT6_KO_B6.2,H3K4me3_E9,H3K4me3_B6.2/all1867_PRMT6_DB_E9_B6.1'
                   '/PRMT6_KO_E.9,PRMT6_KO_B6.2,H3K4me3_E9,H3K4me3_B6.2all.txt', header=0, sep='\t')
#density_based_motif_comparision(peak_df, 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered')
print len(peak_df)
peak_df = peak_df[peak_df['cluster'].isin([0,1,2,4,7,8])]
print len(peak_df)
sample_dict = {'PRMT6_K4me3_occupied': peak_df}
seq = seqOperations.seq4motif(sample_dict)
db = ["JASPAR_CORE_2016_vertebrates.meme", "HOCOMOCOv9.meme", "SwissRegulon_human_and_mouse.meme"]
seqOperations.motif_analysis(db, 10, seq)
'''

### HOMER motif analysis
'''
outpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/HOMERmotif'
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/H3R2_PRMT6_E9_B6_diff.txt', sep='\t', header=0)
seqOperations.motif_analysis_homer(peak_df, path=outpath, samplename='H3R2_PRMT6_E9_B6_diff')
'''

### Performing motif and CpG analysis on prmt6 sites wrt regions eg. TSS, Exon
'''
#df1 = read_csv('/ps/imt/e/HL60_Christene/csv/A1-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered.csv', sep='\t')
#df2 = read_csv('/ps/imt/e/HL60_Christene/further_analysis/overlap/A1-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered_vs_CW4-HL60-rabbit-anti-Ski vs A2-HL60-rabbitIgG filtered_vs_B4-HL60-anti-Runx vs B5-HL60-IgG_matching_Runx filtered.txt', sep='\t')
#uniqueDF = cal_genomic_region.non_overlapping_peaks(df1, df2)

sample_dict = {}
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2ame2_E9 vs IgG_E.9 filtered/H3R2ame2_E9 vs IgG_E.9 filtered.txt',
    header=0, sep='\t')
sample_dict['H3R2me2_E9_50bp'] = peak_df
seq = seqOperations.seq4motif(sample_dict, seqLength='small')
db = ["JASPAR_CORE_2016_vertebrates.meme", "HOCOMOCOv9.meme", "SwissRegulon_human_and_mouse.meme"]
seqOperations.motif_analysis(db, 10, seq)
'''

### Annotate peaks with 6 nearest genes (+,-) strand
'''
diffpeaks = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/PRMT6_new+old_peaks/PRMT6_old+new_only_Enhancer_bound.txt',
    header=0, sep='\t')
gtf_path = '/ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf'
next5genes_annotator(diffpeaks, gtf_path)
'''

### calculate overlaps between peaks
'''
overlapping_samples = {}

sample_name1 = ['H3R2ame2_E9 vs IgG_E.9 filtered']#, 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered']
sample_name2 = ['PRMT6_E9_commer vs IgG_E.9 filtered']#, 'Sample_18F3_RA vs IgG_RA_seq6 filtered']

if len(sample_name1) != len(sample_name2):
    raise ValueError("Unequal sample list for comparison")
else:
    for i in range(0, len(sample_name1)):
        overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, sample_name1[i], sample_name2[i])
        name = overlapping_res.keys()[0]
        with open(basepath + "/further_analysis/overlap/overlapping_peaks.txt", "a") as file:
            file.write(
                name.split('vs')[0] + '\t' + name.split('vs')[2][1:] + '\t' + str(len(overlapping_res.get(name))) + '\n')
'''

## Calculate overlap for external samples
'''
## Put files needed for analysis into 4overlap_external folder before running the anaylsis
sample_name = ['H3R2ame2_E9 vs IgG_E.9 filtered_vs_H3R2me2a_E9_RA vs IgG_E9_RA filtered', 'Sample_18F3 vs Sample_8C9 filtered_vs_Sample_18F3_RA vs IgG_seq6_RA filtered']

filtered_peak_data = {}
for a in sample_name:
    df = read_csv(
        basepath+'/4ovarlap_external/' + a + '.txt', header=0, sep='\t')
    df = df.rename(columns={'Next Gene name': 'Next transcript gene name'})
    filtered_peak_data[a] = df

peakAnalysis_df_ext = {}
for k, v in filtered_peak_data.iteritems():
    name = k
    df = v
    GR_analysis = cal_genomic_region.PeaksAnalysis(df, name, basepath+'/4ovarlap_external/')
    GR_analysis.genomic_regions()
    peakAnalysis_df_ext[name] = GR_analysis
    #GR_analysis.plot_factors('Next Gene biotype')
    #cal_genomic_region.stacke_plot_multiple([name], filtered_peak_data)

sample_name1 = ['H3R2ame2_E9 vs IgG_E.9 filtered_vs_H3R2me2a_E9_RA vs IgG_E9_RA filtered']
sample_name2 = ['Sample_18F3 vs Sample_8C9 filtered_vs_Sample_18F3_RA vs IgG_seq6_RA filtered']
if len(sample_name1) != len(sample_name2):
    raise ValueError("Unequal sample list for comparison")
else:
    for i in range(0, len(sample_name1)):
        overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df_ext, sample_name1[i], sample_name2[i])
        name = overlapping_res.keys()[0]
        with open(basepath+"/further_analysis/overlap/overlapping_peaks.txt", "a") as file:
            file.write(
                sample_name1[i] + '\t' + sample_name2[i] + '\t' + str(len(overlapping_res.get(name))) + '\n')
'''

### Enhancer RNA analysis
'''
from rna_seq.some_methods import enhancer_rna_analysis

peakdf = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/PRMT6_new+old_peaks/PRMT6_old+new_only_Enhancer_bound.txt',
    header=0, sep='\t')
print peakdf.shape
peakdf = peakdf[(peakdf['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'intergenic') | (peakdf['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == 'upstream')]
print peakdf.shape
path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/alignedLane/'
bamRNA = {'NT2D1_E9_1':path+'NT2D1_E9_1/NT2D1_E9_1_GRCh37.bam',
          'NT2D1_E9_2':path+'NT2D1_E9_2/NT2D1_E9_2_GRCh37.bam',
          'NT2D1_E9_3':path+'NT2D1_E9_3/NT2D1_E9_3_GRCh37.bam',
          'NT2D1_B6_1':path+'NT2D1_B6_1/NT2D1_B6_1_GRCh37.bam',
          'NT2D1_B6_2':path+'NT2D1_B6_2/NT2D1_B6_2_GRCh37.bam',
          'NT2D1_B6_3':path+'NT2D1_B6_3/NT2D1_B6_3_GRCh37.bam'}
print(enhancer_rna_analysis(peakdf, bamRNA, outpath=basepath + '/further_analysis/PRMT6_KO_analysis/peak_selecetion_B6/PRMT6_new+old_peaks', samplename='eRNA_PRMT6_n+o_enhancer_bound'))
'''