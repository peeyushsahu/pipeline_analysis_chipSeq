__author__ = 'peeyush'
import subprocess, os
import annotate as Annotate
import pandas as pd


class path():
    def __init__(self):
        basepath = os.path.abspath(os.getcwd())
        #print(basepath)
        self.basepath = os.sep.join(basepath.split(os.sep)[:-2])
        print self.basepath


def create_odir():
    odir = 'results'
    indir = ['alignedLane', 'cache', 'peaks']
    cdpath = path().basepath + '/further_analysis'##os.getcwd() when run script from the folder of interest
    ensure_path(os.path.join(cdpath, odir))
    for i in indir:
        Path = os.path.join(cdpath, odir, i)
        ensure_path(Path)
    print os.path.join(cdpath, odir)
    return os.path.join(cdpath, odir)


def ensure_path(path):
    if not os.path.exists(path):
        os.makedirs(path)


def create_gtf2transcriptDB(gtf_file_path, path, feature='transcript'):
    '''
    This Function will create a database which contains transcript names and their start and end position on chromosome.
    :param gtffile: path to gtf file.
    :return:
    '''
    ## Preparing gtf file for parsing for gene names
    import re
    gtf_file = Annotate.parse_gtf(gtf_file_path)
    gtf_group = gtf_file.groupby("feature")
    gene_gtf = gtf_group.get_group(feature)
    column = ['pid', 'gene_name', 'gene_id', 'gene_biotype', 'transcript_id']
    #gene_gtf.drop('further_information')
    with open(os.path.join(path,'transcript_annotation.db'), 'a') as annotationDB:
        annotationDB.write('chr'+'\t'+'start'+'\t'+'stop'+'\t'+'strand'+'\t'+'gene_name'+'\t'+'pid'+'\t'+'gene_id'+'\t'+'gene_biotype'+'\t'+'transcript_id\n')
        for ind, row in gene_gtf.iterrows():
            gene = '-'
            gene_id = '-'
            gene_biotype = '-'
            pid = '-'
            transcript_id = '-'
            f_cols = row['further_information'].split(";")
            for col in f_cols:
                for name in column:
                    if name in col:
                        if name == 'gene_name': gene = col.split('"')[-2]
                        if name == 'gene_id': gene_id = col.split('"')[-2]
                        if name == 'pid': pid = col.split('"')[-2]
                        if name == 'gene_biotype': gene_biotype = col.split('"')[-2]
                        if name == 'transcript_id': transcript_id = col.split('"')[-2]
            annotationDB.write(str(row['chr'])+'\t'+str(row['start'])+'\t'+str(row['stop'])+'\t'+row['strand']+'\t'+gene+'\t'+pid+'\t'+gene_id+'\t'+gene_biotype+'\t'+transcript_id+'\n')
    trdb_path= os.path.join(path,'transcript_annotation.db')
    create_longestTranscriptDB(trdb_path) # Def for longest transcript
    #########################################


