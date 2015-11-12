__author__ = 'sahu'
import HTSeq.scripts.count as htseq
import pandas as pd
import subprocess as sp

class experiment():
    def __init__(self, name, paths, count_data):
        self.name = name
        self.bampaths = paths
        self.count_data = count_data


def count_data(bam_files, gtf_file, samtype='bam', order='name', stranded='yes', minaqual=10, feature_type='exon', id_attribute='gene_name', overlap_mode='union', samout="", quiet=True):
    '''
    To count exon enrichment for bam using HTSeq-count
    :return:
    '''
    counts_dict = {}
    for bam_file in bam_files:
        print bam_file
        count_dict = htseq.count_reads_in_features(bam_file, gtf_file, samtype, order, stranded, overlap_mode, feature_type, id_attribute, quiet, minaqual, samout)
        name = bam_file.split('/')[-1]
        name = name.split('__')[0]
        counts_dict[name] = count_dict
    count_df = pd.DataFrame(counts_dict)
    count_df_fin = count_df[count_df.sum(axis=1) > 1]
    file = open('/ps/imt/e/HL60_Christene/further_analysis/RNA_seq/GeneId_HL60_SKI_EGFP_KO_stats.txt', 'w')
    file.write('Original feature counts: '+str(len(count_df)))
    file.write('Filtered feature counts: '+str(len(count_df_fin)))
    count_df_fin.to_csv('/ps/imt/e/HL60_Christene/further_analysis/RNA_seq/Filtered_GeneId_HL60_SKI_EGFP_KO.txt', sep='\t')
    #return pd.DataFrame(counts_dict), counts_dict

'''
gtf_file = '/ps/imt/f/Genomes/Homo_sapiens.GRCh37.75.gtf'
bam_files = ['/ps/imt/e/HL60_Christene/results/AlignedLane/HL60_2_10_1__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37/aligned_unique_HL60_2_10_1__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37.bam',
              '/ps/imt/e/HL60_Christene/results/AlignedLane/HL60_2_10_2__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37/aligned_unique_HL60_2_10_2__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37.bam',
              '/ps/imt/e/HL60_Christene/results/AlignedLane/HL60_2_10_3__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37/aligned_unique_HL60_2_10_3__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37.bam',
              '/ps/imt/e/HL60_Christene/results/AlignedLane/HL60_GFP3_1__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37/aligned_unique_HL60_GFP3_1__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37.bam',
              '/ps/imt/e/HL60_Christene/results/AlignedLane/HL60_GFP3_2__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37/aligned_unique_HL60_GFP3_2__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37.bam',
              '/ps/imt/e/HL60_Christene/results/AlignedLane/HL60_GFP3_3__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37/aligned_unique_HL60_GFP3_3__aligned_with_STAR_against_EnsemblGenome_Homo_sapiens_74_37.bam']
'''


def cuffDiff(alignmentLanes, groupA, groupB, genome, label):
    '''
    This will perform cuffdiff of aligned bam files (must be splice aware aligned)
    :param alignmentLanes:
    :param groupA: list of sample names in group A
    :param groupB: list of sample names in group B
    :param genome: Pass the genome object
    :param label: Names of groups e.g. label = ['condition A', 'condition B']
    :return:
    '''
    program = '/home/sahu/Documents/aligners/cufflinks-2.2.1.Linux_x86_64/cuffdiff'
    control = []
    condition = []
    for i in range(0, len(groupA)):
        control.append(control+alignmentLanes[groupA[i]].sortbampath)
        condition.append(condition+alignmentLanes[groupB[i]].sortbampath)
    thread = '-p 6'
    lables = '-L '+','.join(label)
    frag_bias = '-b '+genome.refgenome+'/genome.fa'
    gtfFile = '/ps/imt/f/Genomes/cufflinks_gtf/cuffcmp.combined.gtf'
    outFile = ''
    cmd = ' '.join([program, thread, lables, frag_bias, outFile, gtfFile, control, condition])



def run_cuffDiff(cmd):
    '''
    Runs CuffDiff
    :param cmd:
    :return:
    '''
    print 'Running cuffdiff for '+cmd
    try:
         proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
         stdout, stdrr = proc.communicate()
         print stdrr
         proc.wait()
    except:
        raise IOError ('Subprocess Tophat2 exited with error:', proc)