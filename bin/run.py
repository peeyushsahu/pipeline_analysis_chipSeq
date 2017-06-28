
import alignment.lanes
import alignment.aligner
import rna_seq.bam_processing as bamPrcessing

genome = alignment.aligner.human_GRCh37_74()


raw_lanes = [
    #alignment.lanes.Lane('NT2D1_E9_1', '/ps/imt/f/20151127_RNA/Sample_E_9_C1_091115_R13'),
    #alignment.lanes.Lane('NT2D1_E9_2', '/ps/imt/f/20151127_RNA/Sample_E_9_C2_250915_R8'),
    #alignment.lanes.Lane('NT2D1_E9_3', '/ps/imt/f/20151127_RNA/Sample_E_9_C3_250915_R9'),
    #alignment.lanes.Lane('NT2D1_B6_1', '/ps/imt/f/20151127_RNA/Sample_B6_2_C1_091115_R14'),
    #alignment.lanes.Lane('NT2D1_B6_2', '/ps/imt/f/20151127_RNA/Sample_B6_2_C2_091115_R15'),
    #alignment.lanes.Lane('NT2D1_B6_3', '/ps/imt/f/20151127_RNA/Sample_B6_2_C3_091115_R16'),
    #alignment.lanes.Lane('NT2D1_E9_1_RA', '/ps/imt/f/20160128_RNA/Sample_R17_E9_1p_141215'),
    #alignment.lanes.Lane('NT2D1_E9_2_RA', '/ps/imt/f/20160128_RNA/Sample_R18_E9_2p_141215'),
    #alignment.lanes.Lane('NT2D1_E9_3_RA', '/ps/imt/f/20160128_RNA/Sample_R19_E9_3p_141215'),
    #alignment.lanes.Lane('NT2D1_B6_1_RA', '/ps/imt/f/20160128_RNA/Sample_R20_B62_1p_141215'),
    #alignment.lanes.Lane('NT2D1_B6_2_RA', '/ps/imt/f/20160128_RNA/Sample_R21_B62_2p_141215'),
    #alignment.lanes.Lane('NT2D1_B6_3_RA', '/ps/imt/f/20160128_RNA/Sample_R22_B62_3p_141215'),

    alignment.lanes.Lane('Histone_3_RA', '/ps/imt/f/20140114/Sample_H3+'),
]

# Aliging read files with chosen aligner
# Aligner = Bowtie2, Tophat2
alignedLane = {}
bamPaths = []
for lanes in raw_lanes:
    lanes.join_multiple_fq()
    lanes.do_alignment(genome, 'Bowtie2')
    bamPaths.append(lanes.bampath)

# generating count for features from bam files

#bamPrcessing.count_Data_featureCounts(bamPaths, genome.gtfFile, count_out='/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_+ATRA/count/NT2D1_3d_+ATRA_tophat2_KO.txt')

#  RNASeq diffcalling CuffDiff
#controlSamples = ['NT2D1_E9_1_RA', 'NT2D1_E9_2_RA', 'NT2D1_E9_3_RA']
#conditionSamples = ['NT2D1_B6_1_RA', 'NT2D1_B6_2_RA','NT2D1_B6_3_RA']
#bamPrcessing.cuffDiff(alignedLane, controlSamples, conditionSamples, genome, ['PRMT6_KO_EGFP9_RA', 'PRMT6_KO_B6_RA'])

