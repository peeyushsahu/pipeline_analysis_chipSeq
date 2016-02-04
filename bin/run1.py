import os

from pandas import read_csv
from alignment import retroT_analysis
from annotate.Annotate import next5genes_annotator

from overlap_analysis import cal_genomic_region, differential_binding, filterPeaks
from overlap_analysis.differential_binding import modification4nearestgenes
from plotsAndseq import seqOperations
import plotsAndseq.plots as plots
from plotsAndseq.plots import GR_heatmaps_DF_for_peaks
from plotsAndseq.seqOperations import density_based_motif_comparision
import pandas as pd
#import pdb; pdb.set_trace()
import gc


__author__ = 'peeyush'


#time.strftime("%d/%m/%y")

folders = ["overlap",
           "differential",
           "filtered",
           "plots",
           "seq4motif",
           "valley_peaks",
           "CpG",
           "density_based_motif"]
for folder in folders:
    #print 'Directory_for_result: ' + '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/'+folder
    path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/'+folder
    if not os.path.exists(path):
        os.makedirs(path)
print 'Output folder created'

sample_name = [
     #'H3K4me3_B5.1 vs IgG_B5.1 filtered',
     #'Sample_K27ac vs IgG_seq6 filtered',
     #'Sample_K27ac_10_wt vs IgG_seq6_RA filtered',
     #'PRMT6_seq5 vs IgG_seq4 filtered',
     #'Sample_K4me1_RA vs Sample_PIS_RA filtered',
     #'Sample_PRMT6_3_RA vs IgG_seq6_RA filtered',
     #'H3K4me3_seq2 vs IgG_seq2 filtered',
     #'Sample_EZH1 vs IgG_seq6 filtered',
     #'PRMT6_seq6_RA vs IgG_seq6_RA filtered',
     #'H3K4me3_seq2 vs Sample_PIS filtered',
     #'Sample_EZH2_RA vs IgG_seq6_RA filtered',
     #'PRMT6_KO_B5.1 vs IgG_B5.1 filtered',
     #'Sample_K27ac_RA vs IgG_seq6_RA filtered',
     #'Sample_K9me3_RA vs IgG_seq6_RA filtered',
     #'PRMT6_KO_10.8 vs IgG_10.8 filtered',
     #'H3K4me3_E9 vs IgG_E.9 filtered',
     #'Sample_18F3 vs Sample_8C9 filtered',
     #'Sample_K9me3 vs IgG_seq6 filtered',
     #'Sample_pol2_RA vs IgG_seq6_RA filtered',
     #'PRMT6_3_seq1 vs IgG_seq2 filtered',
     #'Sample_PRMT6_2_10_wt vs IgG_seq6_RA filtered',
     #'Sample_K27me1_RA vs IgG_seq6_RA filtered',
     #'YY1_seq3 vs IgG_seq4 filtered',
     #'Sample_K36me3_RA vs IgG_seq6_RA filtered',
     #'Sample_EZH2 vs IgG_seq6 filtered',
     #'H3K27ac_B5.1 vs IgG_B5.1 filtered',
     #'H3K4me3_seq2_RA vs Sample_PIS_RA filtered',
     #'H3K27ac_E9 vs IgG_E.9 filtered',
     #'PRMT6_seq6 vs IgG_seq6 filtered',
     #'Sample_PRMT6_3_0_wt vs IgG_seq6 filtered',
     #'Sample_K27me1 vs IgG_seq6 filtered',
     #'YY1_seq3_RA vs IgG_seq4_RA filtered',
     #'Sample_K4me1 vs Sample_PIS filtered',
     #'Sample_EZH1_RA vs IgG_seq6_RA filtered',
     #'PRMT6_1_seq1 vs IgG_seq2 filtered',
     #'Sample_K27me3_RA vs IgG_seq6_RA filtered',
     #'H3R2ame2_E9 vs IgG_E.9 filtered',
     #'Sample_K27me3 vs IgG_seq6 filtered',
     #'Sample_K36me3 vs IgG_seq6 filtered',
     #'Sample_pol2 vs IgG_seq6 filtered',
     #'Sample_K4me3_10_wt vs IgG_seq6_RA filtered',
     #'Sample_18F3_RA vs IgG_seq6_RA filtered',
     #'PRMT6_KO_B6.2 vs IgG_B6.2 filtered',
     #'H3R2ame2_B5.1 vs IgG_B5.1 filtered',
     #'PRMT6_seq5_RA vs IgG_seq4_RA filtered',
     #'PRMT6_KO_E.9 vs IgG_E.9 filtered',
     #'PRMT6_seq4 vs IgG_seq4 filtered'
     ]


# Here import peak called data in a list....
Filtering = True
if Filtering:
    peak_data = {}
    for a in sample_name:
        df = read_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/' + a + '.csv',
            header=0, sep='\t')
        df = df.rename(columns={'Next Gene name':'Next transcript gene name'})
        peak_data[a] = df
    print "Number of sample are being analysed: ", peak_data.__len__()

    print "Filtering peaks."
    filtered_peak_data = filterPeaks.filterpeaks(peak_data)
else:
    filtered_peak_data = {}
    print "############ No filtering ###############"
    for a in sample_name:
        df = read_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/csv/' + a + '.csv',
            header=0, sep='\t')
        df = df.rename(columns={'Next Gene name': 'Next transcript gene name'})
        filtered_peak_data[a] = df
    print "Number of sample are being analysed: ", filtered_peak_data.__len__()

'''
df = read_csv(
            '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/PRMT6_DB_E9_B5.1.csv',
            header=0, sep='\t')
filtered_peak_data['PRMT6_peaks'] = df
## Plot stacked plot for selected samples
cal_genomic_region.stacke_plot_multiple(['PRMT6_peaks']
                                        , filtered_peak_data)
cal_genomic_region.stacke_plot_multiple(['H3K4me3_seq2 vs IgG_seq2 filtered', 'H3K4me3_RA_seq2 vs IgG_RA_seq2 filtered']
                                        , filtered_peak_data)
cal_genomic_region.stacke_plot_multiple(['Sample_18F3 vs Sample_8C9 filtered', 'Sample_18F3_RA vs IgG_RA_seq6 filtered']
                                        , filtered_peak_data)
'''
peakAnalysis_df = {}
for k, v in filtered_peak_data.iteritems():
    name = k
    df = v
    GR_analysis = cal_genomic_region.PeaksAnalysis(df, name)
    GR_analysis.genomic_regions()
    peakAnalysis_df[name] = GR_analysis
    GR_analysis.plot_factors('Next Gene biotype')


### Performs differential binding calulation from full sample
'''
sample = ['PRMT6_seq6', 'PRMT6_KO_B6.2', 'PRMT6_KO_B5.1']
diffbind = differential_binding.Overlaps(sample, filtered_peak_data)
diffbind.diffBinding('PRMT6_seq6 vs IgG_seq6 filtered')
'''
### Diff.binding for nearest genes
'''
sample = ['Sample_K36me3', 'Sample_K36me3_RA',
          'Sample_pol-2', 'Sample_pol-2_RA',
          'H3K4me3_seq2', 'H3K4me3_RA_seq2',
          'Sample_K27me3', 'Sample_K27me3_RA']
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/diffPeaks_P6_nearest_5_genes.csv',
    header=0, sep=',')

modification4nearestgenes(peak_df, 'prmt6_nearest5genes', sample)
'''

### Diff. Binding calculation from altered sample (external)

#'H3R2ame2_E9', 'H3K36me3_E9', 'H3K36me3_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2', 'H3K27ac_E9','H3K27ac_B6.2'
multiple_df = ['/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/totalDEgenes_analysis/chip_H3K4me3_down.csv',
               '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/totalDEgenes_analysis/chip_H3K4me3_up.csv',
               '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/totalDEgenes_analysis/chip_H3K27ac_down.csv',
               '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/totalDEgenes_analysis/chip_H3K27ac_up.csv',
               '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/totalDEgenes_analysis/chip_H3K4me1_up.csv',
               '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/totalDEgenes_analysis/chip_H3K4me1_down.csv']
for df in multiple_df:
    sample = ['H3K4me3_E9', 'H3K4me3_B6.2', 'H3K27ac_E9','H3K27ac_B6.2', 'H3K4me1_E9', 'H3K4me1_B6']
    peak_df = read_csv(df, header=0, sep='\t')
    filtered_peak = {'loaded_sample': peak_df}
    diffbind = differential_binding.Overlaps(sample, filtered_peak)
    diffbind.diffBinding('loaded_sample', outpath=df+'_diff.txt') #, genewide=True

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
### Gene-wide chip profile for broad histone marks
'''
bam_list = ['H3K36me3_E9', 'H3K36me3_B6.2']
sample_name = 'R2_K36_all946'
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/H3R2me2a/H3R2ame2_E9,H3K36me3_E9,H3K36me3_B6.2,H3K4me3_E9,H3K4me3_B6.2,H3K27ac_E9,H3K27ac_B6.2'
                   '/diff_binding_all946.txt', header=0, sep='\t')
#peak_df = peak_df[peak_df['log2FC_H3K27ac_E9_norm_vs_H3K27ac_B6.2_norm'] > 2]
plots.grHeatmap4wholeGene(peak_df, bam_list, sample_name)
'''

### Compares multiple ChIP-Seq profile using peaks (heatmap) from on sample
'''
region = ['all'] #'all', 'tss', 'exon', 'intron', 'intergenic', 'upstream'
bam_list = [['PRMT6_seq4', 'PRMT6_seq6', 'PRMT6_KO_E.9', 'PRMT6_KO_B6.2']]
for bams in bam_list:
    for i in region:
        GR_heatmaps_DF_for_peaks(bams, filtered_peak_data.get('PRMT6_seq4 vs IgG_seq4 filtered'), region=i,
                                 sort=False, sort_column='H3R2ame2_E9', scale_df=False)
'''
### Comapre ChIP-Seq profile from altered sample (external)
'''
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/filtered/'
                   'H3R2ame2_E9 vs IgG_E.9 filtered.csv', header=0, sep=',')
print peak_df.head(2)
## filter peak data with matching column using list
listGroups = ['unchangedPeaks','upRegulatedPeaks','downRegulatedPeaks']
for List in listGroups:
    List = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/'
                    'H3R2me2a/'+List+'.txt', sep='\t', header=None)
    peaks_df_n = filterPeaks.extract_data(peak_df, List.loc[:,0].tolist(), 'Next transcript gene name')
    print peaks_df_n.head()
    #peaks_df_n.to_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/multidimensional_analysis/'
    #                  'H3R2me2a/'+List+'peakdf.txt', sep='\t')
    # bam_list = ['PRMT6_seq6','PRMT6_seq6_RA','Sample_K27ac','Sample_K27ac_RA','Sample_K4me1','Sample_K4me1_RA',
    # 'PRMT6_KO_E.9','PRMT6_KO_B5.1','H3K27ac_E9','H3K27ac_B5.1',
    # 'H3R2ame2_E9', 'H3R2me2a_B6.2', 'H3K4me3_E9','H3K4me3_B6.2']
    bam_list = [['H3R2ame2_E9', 'H3K36me3_E9', 'H3K36me3_B6.2', 'H3K4me3_E9', 'H3K4me3_B6.2', 'H3K27ac_E9',
                 'H3K27ac_B6.2']]

    ### If DF is from R change column names ('' ','.')
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
            #peak_df['cluster'] = d.Categorical(peak_df['cluster'], [6,7,4,3,0,1,5,8,2]) #7,5,4,0,1,2,3,9,8,6
            #peak_df = peak_df.sort('cluster')
            #peak_df = peak_df[(peak_df['cluster'] == 4) | (peak_df['cluster'] == 5) | (peak_df['cluster'] == 7)]
            GR_heatmaps_DF_for_peaks(List, peaks_df_n, region=i, sort=False, sort_column='H3K4me3_seq2', scale_df=False)
            #GPcount = peak_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'].value_counts()
            #GPcount = zip(GPcount.index, GPcount.values)
            #cal_genomic_region.plotGenomicregions(GPcount, 'DiffBind_P6_vs_all')
    gc.collect()
'''
### Density based motif analysis
'''
peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlapping_plots/PRMT6_2_seq6,Sample_pol-2,Sample_K36me3,H3K4me3_seq2,Sample_K27me1,Sample_K27ac,Sample_K4me1,Sample_K9me3,Sample_K27me3,Sample_18F3/intergenic7325'
                   '/PRMT6_2_seq6,Sample_pol-2,Sample_K36me3,H3K4me3_seq2,Sample_K27me1,Sample_K27ac,Sample_K4me1,Sample_K9me3,Sample_K27me3,Sample_18F3intergenic.csv',
    header=0, sep=',')
#density_based_motif_comparision(peak_df, 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered')
peak_df = peak_df[peak_df['cluster'] == 7]
sample_dict = {'PRMT6_intergenic_vs_H3K9_Non_RA': peak_df}
seq = seqOperations.seq4motif(sample_dict)
db = ["JASPAR_CORE_2014_vertebrates.meme", "uniprobe_mouse.meme"]
seqOperations.motif_analysis(db, 10, seq)
'''
### Performing motif and CpG analysis on prmt6 sites wrt regions eg. TSS, Exon
'''
sample_dict = {}
prmt6_df = filtered_peak_data.get('PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered')
for i in ['tss', 'exon', 'intron', 'intergenic', 'upstream']:
    sample_dict['PRMT6_2_RA_seq6_'+i] = prmt6_df[prmt6_df['GenomicPosition TSS=1250 bp, upstream=5000 bp'] == i]
seq = seqOperations.seq4motif(sample_dict)
db = ["JASPAR_CORE_2014_vertebrates.meme", "uniprobe_mouse.meme"]
seqOperations.motif_analysis(db, 10, seq)
'''
'''
sample_dict = {}
prmt6_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/H3R2ame2_E9_Sample_18F3_H3R2ame2_B5.1.csv',
    header=0, sep=',')
sample_dict['H3R2me2a_B5.1_KO'] = prmt6_df
seq = seqOperations.seq4motif(sample_dict)
db = ["JASPAR_CORE_2016_vertebrates.meme", "HOCOMOCOv9.meme", "SwissRegulon_human_and_mouse.meme"]
seqOperations.motif_analysis(db, 10, seq)
'''
### Annotate peaks with 5 nearest genes (+,-) strand
'''
diffpeaks = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential/DiffBind_P6_vs_all.csv',
    header=0, sep=',')
gtf_path = '/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf'
next5genes_annotator(diffpeaks, gtf_path)
'''
### calculate overlaps between peaks
'''
overlapping_samples = {}

sample_name1 = ['H3K4me3_seq2 vs IgG_seq2 filtered']#, 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered']
sample_name2 = ['H3K4me3_E9 vs IgG_E.9 filtered']#, 'Sample_18F3_RA vs IgG_RA_seq6 filtered']

if len(sample_name1) != len(sample_name2):
    raise ValueError("Unequal sample list for comparison")
else:
    for i in range(0, len(sample_name1)):
        overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, sample_name1[i], sample_name2[i])
        name = overlapping_res.keys()[0]
        with open("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/overlapping_peaks.txt", "a") as file:
            file.write(
                name.split('vs')[0] + '\t' + name.split('vs')[2][1:] + '\t' + str(len(overlapping_res.get(name))) + '\n')
'''

## caluculate peaks overlap using external data
'''
overlapping_samples = {}
sample_name3 = ['PRMT6_DB_E9_B.1']
sample_name4 = ['filteredYY1_seq3 vs IgG_seq2 filtered']

diff_sample = ['PRMT6_DB_E9_B5.1', 'filteredYY1_seq3 vs IgG_seq2 filtered'] ## load samples to analysis

#peak_data2 = {}
for a in diff_sample:
    df = read_csv(
        '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/' + a + '.csv',
        header=0, sep='\t')
    names_list = a.split('/')
    if len(names_list) > 1:
        a = names_list[len(names_list)-1]
    GPcount = df['GenomicPosition TSS=1250 bp, upstream=5000 bp'].value_counts()
    GPcount = zip(GPcount.index, GPcount.values)
    cal_genomic_region.plotGenomicregions(GPcount, a)
    peakAnaObj = cal_genomic_region.PeaksAnalysis(peaks_df=df, con_name=a)
    overlapping_samples[a] = peakAnaObj
    print a

for i in range(0, len(sample_name3)):
    print 'External sample overlapping:', sample_name3[i], '======', sample_name4[i]
    overlap4diff = cal_genomic_region.OverlappingPeaks(overlapping_samples, sample_name3[i], sample_name4[i])
    name = overlap4diff.keys()[0]
    with open("/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/overlapping_peaks.txt", "a") as file:
        file.write(
            name + '\t' + str(len(overlap4diff.get(name))) + '\n')
    overlapping_samples[overlap4diff.keys()[0]] = pd.DataFrame(overlap4diff.get(overlap4diff.keys()[0]))
'''
#this will make seq from Genomic regions
#print "Samples for seq fetching:", len(overlapping_samples)
#seq = seqOperations.seq4motif(filtered_peak_data)
#db = ["JASPAR_CORE_2014_vertebrates.meme", "uniprobe_mouse.meme"]
#seqOperations.motif_analysis(db, 10, seq)

#seq1 = seqOperations.seq4motif({'p6_H3K9':peak_df})
#db1 = ["JASPAR_CORE_2014_vertebrates.meme", "uniprobe_mouse.meme"]
#seqOperations.motif_analysis(db1, 10, seq1)
