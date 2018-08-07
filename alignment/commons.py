__author__ = 'peeyush'
import subprocess as sp
import os
import annotate as Annotate
import pandas as pd


class path():
    def __init__(self):
        basepath = os.path.abspath(os.getcwd())
        #print(basepath)
        self.basepath = os.sep.join(basepath.split(os.sep)[:-2])
        print(self.basepath)


def create_odir():
    odir = 'results'
    indir = ['alignedLane', 'cache', 'peaks']
    cdpath = path().basepath + '/further_analysis'##os.getcwd() when run script from the folder of interest
    ensure_path(os.path.join(cdpath, odir))
    for i in indir:
        Path = os.path.join(cdpath, odir, i)
        ensure_path(Path)
    print(os.path.join(cdpath, odir))
    return os.path.join(cdpath, odir)


def ensure_path(path):
    if not os.path.exists(path):
        os.makedirs(path)


def to_bed(peaks):
    bed = []
    for ind, row in peaks.iterrows():
        bed.append(['chr'+str(row['chr']), int(row['start']), int(row['stop'])])
    return pd.DataFrame(bed)


def peakdf_columns():
    '''
    Minimum column requirement for analysis
    :return:
    '''
    columns_2_select = [
        'chr',
        'start',
        'stop',
        'GenomicPosition TSS=1250 bp, upstream=5000 bp',
        'Next Transcript tss distance',
        'Next transcript gene name',
        'Next Transcript stable_id',
        'Next transcript strand',
        'summit'
    ]
    return columns_2_select

#'NT2D1_KO_H3R2me2a_-Doxy_R1': 'aligned_unique_H2R2me2a_B6_rescue__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_H3R2me2a_+Doxy_R1': 'aligned_unique_H2R2me2a_B6_rescue_Dox__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_H3R2me2a_-RA_R1': 'aligned_unique_H3R2me2a_r2_B6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_H3R2me2a_-Doxy_R2': 'aligned_unique_H3R2me2a_r2_B6_rescue__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_H3R2me2a_+Doxy_R2': 'aligned_unique_H3R2me2a_r2_B6_rescue_Dox__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3R2me2a_-RA_R1': 'aligned_unique_H3R2me2a_r2_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3R2me2a_-Inh_R1': 'aligned_unique_H3R2me2a_r3_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3R2me2a_+Inh_R1': 'aligned_unique_H3R2me2a_r3_E9_inh_16ug__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_H3R2me2a_-RA_R2': 'aligned_unique_H3R2me2a_r4_B6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3R2me2a_-RA_R2': 'aligned_unique_H3R2me2a_r4_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3R2me2a_-Inh_R2': 'aligned_unique_H3R2me2a_r4_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3R2me2a_+Inh_R2': 'aligned_unique_H3R2me2a_r4_E9_inh_16ug__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_KMT2D_-RA_abgent': 'aligned_unique_MLL4_abgen_B6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_KMT2D_+RA_abgent': 'aligned_unique_MLL4_abgen_B6_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_KMT2D_-RA_abgent': 'aligned_unique_MLL4_abgen_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_KMT2D_+RA_abgent': 'aligned_unique_MLL4_abgen_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'U2OS_CT_H3R2me2a': 'aligned_unique_U2OS_H3R2me2a_CT1__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'U2OS_KO_H3R2me2a': 'aligned_unique_U2OS_H3R2me2a_KO2__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_KMT2D_-RA_santacruz': 'aligned_unique_MLL4_santac_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_KMT2D_+RA_santacruz': 'aligned_unique_MLL4_santac_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_KMT2D_-RA_santacruz': 'aligned_unique_MLL4_santac_B6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_KMT2D_+RA_santacruz': 'aligned_unique_MLL4_santac_B6_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3K27me3_-RA': 'aligned_unique_H3K27me3_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_H3K27me3_-RA': 'aligned_unique_H3K27me3_B6__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_CT_H3K27me3_+RA': 'aligned_unique_H3K27me3_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
#'NT2D1_KO_H3K27me3_+RA': 'aligned_unique_H3K27me3_B6_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam',
bam_files = {
    'NT2D1_CT_H3R2me2a_-RA_R2': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3R2me2a_r4_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup/aligned_unique_H3R2me2a_r4_E9__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_dedup.bam'
}


def bam_to_bigwig(bamfiles_dict, outpath=None):
    '''
    Convert bam files to bigwig using deeptools.bamCoverage
    :return:
    '''

    def run_cmd(cmd):
        '''
        Using subprocees
        :param cmd:
        :return:
        '''
        try:
            proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            stdout, stdrr = proc.communicate()
            print(stdrr)
            print(stdout)
        except Exception as e:
            raise IOError('Subprocessexited with error:', e)

    for name, path in bamfiles_dict.items():
        cmd = ['python3']
        cmd.extend(['/ps/imt/tools/deepTools-2.5.4/bin/bamCoverage'])
        cmd.extend(['-b', path])
        cmd.extend(['-of bedgraph'])  # bedgraph or bigwig
        cmd.extend(['--binSize 50'])
        if 'santacruz' in name:
            print('Here')
            cmd.extend(['-o', os.path.join(outpath, name+'normCPM.bg')])
            cmd.extend(['--normalizeUsing CPM'])
            cmd.extend(['--effectiveGenomeSize 2864785220'])
        else:
            cmd.extend(['-o', os.path.join(outpath, name+'.bg')])
        #print(cmd)
        cmd = ' '.join(cmd)
        print(cmd)
        run_cmd(cmd)


def ensembl_bedgraph_to_ucsc(folder_path, outpath=None):
    '''
    Add chr in chromosome names
    :return:
    '''
    for name in os.listdir(folder_path):
        bedgraph_df = pd.read_csv(os.path.join(folder_path, name), header=None, sep='\t')
        print(name)
        print('Size of file:', bedgraph_df.shape)
        bedgraph_df[0] = bedgraph_df[0].astype('str')
        bedgraph_df = bedgraph_df[bedgraph_df[0].str.len() < 4]
        print('Size of file after:', bedgraph_df.shape)
        bedgraph_df = bedgraph_df[~bedgraph_df[0].isin(['MT'])]
        bedgraph_df[0] = 'chr' + bedgraph_df[0]
        bedgraph_df = bedgraph_df.sort_values([0, 1], ascending=[True, True])
        print('Size of file after MT rem:', bedgraph_df.shape)
        print('###################################################')
        bedgraph_df.to_csv(os.path.join(outpath, name), header=False, index=False, sep='\t')


def ucscbedGraph_to_ucscbigwig(folder_path, outpath=None):
    '''
    Convert bam files to bigwig using deeptools.bamCoverage
    :return:
    '''

    def run_cmd1(cmd):
        '''
        Using subprocees
        :param cmd:
        :return:
        '''
        try:
            proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            stdout, stdrr = proc.communicate()
            print(stdrr)
            print(stdout)
            proc.wait()
        except Exception as e:
            raise IOError('Subprocessexited with error:', e)

    for name in os.listdir(folder_path):
        cmd = []
        cmd.extend(['/ps/imt/tools/bedgraph2bigwig/bedGraphToBigWig'])
        cmd.extend([os.path.join(folder_path, name)])
        cmd.extend(['/ps/imt/tools/bedgraph2bigwig/hg19.chrom.sizes'])
        cmd.extend([os.path.join(outpath, str(name).split('.')[0] + '.bw')])
        cmd = ' '.join(cmd)
        print(cmd)
        run_cmd1(cmd)


def prepare_md5sum(folder_path, file_extention=None):
    '''
    Write md5Sum for all files in a given folder
    '''
    folder_path = folder_path
    #for f_name in os.listdir(folder_path):
    with open(os.path.join(folder_path, 'md5sum'), mode='a') as file:
        for name in os.listdir(os.path.join(folder_path)):
            file_path = os.path.join(folder_path, name)
            if os.path.isfile(file_path):
                print(name)
                cmd = ['md5sum', file_path]
                proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
                stdout, stderr = proc.communicate()
                #print(stdout)
                #print(stderr)
                output = stdout.decode('ascii').split(' ')
                file.write(name+'\t'+output[0]+'\n')
                proc.wait()
            else:
                for name1 in os.listdir(file_path):
                    print(name1)
                    sub_file_path = os.path.join(file_path, name1)
                    cmd = ['md5sum', sub_file_path]
                    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
                    stdout, stderr = proc.communicate()
                    #print(stdout)
                    #print(stderr)
                    output = stdout.decode('ascii').split(' ')
                    file.write(name1+'\t'+name+'\t'+output[0]+'\n')
                    proc.wait()
        file.close()

