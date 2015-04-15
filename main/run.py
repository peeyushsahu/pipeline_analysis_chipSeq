from overlap_analysis import cal_genomic_region, filterPeaks

__author__ = 'peeyush'
from pandas import read_csv
import os


paths = ["/ps/imt/e/20150105_AG_Bauer_H3R2/further_analysis", "/ps/imt/e/20150105_AG_Bauer_H3R2/further_analysis/overlap",
        "/ps/imt/e/20150105_AG_Bauer_H3R2/further_analysis/diffrential", "/ps/imt/e/20150105_AG_Bauer_H3R2/further_analysis/filtered",
        "/ps/imt/e/20150105_AG_Bauer_H3R2/further_analysis/plots"]

# create directories
for path in paths:
    if not os.path.exists(path):
        os.makedirs(path)

sample_name = ['Encode_NT2D1_H3K4me1 encode', 'Encode_NT2D1_H3K36me3 encode',
               'PRMT6_2_seq5 vs IgG_seq2 filtered', 'PRMT6_2_seq6 vs IgG_seq6 filtered','PRMT6_2_RA_seq5 vs IgG_RA_seq2 filtered', 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered']

# Here import peak called data in a list....
peak_data = {}
for a in sample_name:
    print a
    df = read_csv('/home/peeyush/Desktop/Analysis_ChipSeq_Genes/Samples/New_B2_MACS_20-OCT-14/All_samples/csv/'+a+'.csv', header=0, sep='\t')
    peak_data[a] = df
print peak_data.__len__()

filtered_peak_data = filterPeaks.filterpeaks(peak_data)

# Then pass this list to genomic regions and annotate
# [name,DF,GR]
peakAnalysis_df = {}
for k, v in filtered_peak_data.iteritems():
    #print type(lists)
    name = k
    df = v
    #print len(df)
    GR_analysis = cal_genomic_region.PeaksAnalysis(df, name)
    cal_genomic_region.genomic_regions(GR_analysis)
    peakAnalysis_df[name] = GR_analysis
#cal_genomic_region.multiPiePlot(peakAnalysis_df)
#print 'its coming here'
sample_name1 = ['PRMT6_2_seq6 vs IgG_seq6 filtered', 'PRMT6_2_seq5 vs IgG_seq2 filtered', 'PRMT6_2_RA_seq5 vs IgG_RA_seq2 filtered', 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered',
                'PRMT6_2_seq6 vs IgG_seq6 filtered', 'PRMT6_2_seq5 vs IgG_seq2 filtered', 'PRMT6_2_RA_seq5 vs IgG_RA_seq2 filtered', 'PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered']
sample_name2 = ['Encode_NT2D1_H3K4me1 encode','Encode_NT2D1_H3K4me1 encode','Encode_NT2D1_H3K4me1 encode', 'Encode_NT2D1_H3K4me1 encode',
                'Encode_NT2D1_H3K36me3 encode','Encode_NT2D1_H3K36me3 encode','Encode_NT2D1_H3K36me3 encode','Encode_NT2D1_H3K36me3 encode']


overlapping_samples = {}
for i in range(0, len(sample_name1)):
    overlapping_res = cal_genomic_region.OverlappingPeaks(peakAnalysis_df, sample_name1[i], sample_name2[i], Encode='Encode')
    #overlapping_obj = overlapping_analysis.Overlaps(overlapping_res, peak_data)
    #overlapping_analysis.diffBinding(overlapping_obj)
    print overlapping_res.keys()[0]
    #print overlapping_res.get(overlapping_res.keys()[0])
    overlapping_samples.update({overlapping_res.keys()[0]: overlapping_res.get(overlapping_res.keys()[0])})

#overlap_list1 = ['PRMT6_2_seq1 vs IgG_seq1 filtered_vs_PRMT6_2_seq3 vs IgG_seq2 filtered']#, 'PRMT6_2_RA_seq2 vs IgG_RA_seq2 filtered_vs_PRMT6_2_RA_seq3 vs IgG_seq2 filtered',
                 #'PRMT6_2_seq5 vs IgG_seq2 filtered_vs_PRMT6_2_seq6 vs IgG_seq6 filtered', 'PRMT6_2_seq5 vs IgG_seq2 filtered_vs_PRMT6_2_seq6 vs IgG_seq6 filtered',
                 #'PRMT6_2_RA_seq2 vs IgG_RA_seq2 filtered_vs_PRMT6_2_RA_seq3 vs IgG_seq2 filtered', 'PRMT6_2_RA_seq2 vs IgG_RA_seq2 filtered_vs_PRMT6_2_RA_seq3 vs IgG_seq2 filtered']

#overlap_list2 = ['PRMT6_2_seq6 vs IgG_seq6 filtered'] #,'PRMT6_2_RA_seq5 vs IgG_RA_seq2 filtered_vs_PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered',
                 #'H3K4me3_seq2 vs IgG_seq2 filtered', 'H3K27me3_seq2 vs IgG_seq2 filtered', 'H3K4me3_RA_seq2 vs IgG_RA_seq2 filtered', 'H3K27me3_RA_seq2 vs IgG_RA_seq2 filtered']
#for i in range(0, len(overlap_list1)):
#    cal_genomic_region.overlappingPeaksLess(overlapping_samples, filtered_peak_data, overlap_list1[i], overlap_list2[i])
#print overlapping_res.get(overlapping_res.keys()[0])
