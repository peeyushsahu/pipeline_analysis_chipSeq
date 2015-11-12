__author__ = 'sahu'
import HTSeq.scripts.count as htseq
import pandas as pd
import subprocess as sp
import alignment.commons as common
import os

class experiment():
    def __init__(self, name, paths, count_data):
        self.name = name
        self.bampaths = paths
        self.count_data = count_data


def count_data(bam_files, gtf_file, resultpath, name = '_', samtype='bam', order='name', stranded='yes', minaqual=10, feature_type='exon', id_attribute='gene_name', overlap_mode='union', samout="", quiet=True):
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
    file = open(os.path.join(resultpath, 'RNA_seq', name, '_htseq_counts_stats.txt'), 'w')
    file.write('Original feature counts: '+str(len(count_df)))
    file.write('Filtered feature counts: '+str(len(count_df_fin)))
    pathcountData = os.path.join(resultpath, 'RNA_seq', name, '_htseq_counts.txt')
    count_df_fin.to_csv(pathcountData, sep='\t')
    return pathcountData


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
    result_dir = None
    program = '/home/sahu/Documents/aligners/cufflinks-2.2.1.Linux_x86_64/cuffdiff'
    control = []
    condition = []
    for i in range(0, len(groupA)):
        control.append(control+alignmentLanes[groupA[i]].sortbampath)
        condition.append(condition+alignmentLanes[groupB[i]].sortbampath)
        if result_dir is None: result_dir = alignmentLanes[groupA[i]].resultdir
    thread = '-p 6'
    lables = '-L '+','.join(label)
    frag_bias = '-b '+genome.refgenome+'/genome.fa'
    gtfFile = '/ps/imt/f/Genomes/cufflinks_gtf/cuffcmp.combined.gtf'
    outFile = '-o '+common.ensure_path(os.path.join(result_dir, 'RNAseq', 'cuffDiff'))
    cmd = ' '.join([program, thread, lables, frag_bias, outFile, gtfFile, control, condition])
    run_cuffDiff(cmd)



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