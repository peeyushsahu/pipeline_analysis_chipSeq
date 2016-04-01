__author__ = 'peeyush'
import subprocess, os
import annotate.Annotate as Annotate
import pandas as pd


class path():
    def __init__(self):
        basepath = os.path.abspath(os.getcwd())
        self.basepath = os.sep.join(basepath.split(os.sep)[:-2])
        print self.basepath


def create_odir():
    odir = 'results'
    indir = ['alignedLane', 'cache', 'peaks']
    cdpath = path().basepath + '/further_analysis'##os.getcwd() when run script from the folder of interest
    if not os.path.exists(os.path.join(cdpath, odir)): os.mkdir(os.path.join(cdpath, odir))
    for i in indir:
        print os.path.join(cdpath, odir, i)
        if not os.path.exists(os.path.join(cdpath, odir, i)):
            os.mkdir(os.path.join(cdpath, odir, i))
    print os.path.join(cdpath, odir)
    return os.path.join(cdpath, odir)


def ensure_path(path):
    if not os.path.exists(path):
        os.mkdir(path)


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

def create_longestTranscriptDB(path):
    '''
    This will filter longest transcript per gene.
    :param path: path of transcript_annotation.db from def create_gtf2transcriptDB()
    :return:
    '''
    ## Parsing gtf file and storing gene names with position.
    gene_gtf = pd.read_csv(os.path.join(path,'transcript_annotation.db'), sep='\t', header=0)
    gene_gtf['chr'] = gene_gtf['chr'].astype('str')
    gene_gtf = gene_gtf[gene_gtf['chr'].str.strip.len() < 4]
    gene_gtf = gene_gtf.sort('gene_name', ascending=True)
    gene_name = None
    chr = ''
    start = None
    stop = None
    strand = ''
    pid = '-'
    gene_id = ''
    gene_biotype = ''
    transcript_id = ''
    length = 0
    with open(os.path.join(path,'longest_transcript_annotation.db'), 'a') as annotationDB:
        annotationDB.write('chr'+'\t'+'start'+'\t'+'stop'+'\t'+'length'+'\t'+'strand'+'\t'+'gene_name'+'\t'+'pid'+'\t'+'gene_id'+'\t'+'gene_biotype'+'\t'+'transcript_id\n')
        for index, rows in gene_gtf.iterrows():
            if gene_name == rows['gene_name']:
                if length < (int(rows['stop']) - int(rows['start'])):
                    #print length
                    length = rows['stop'] - rows['start']
                    gene_name = rows['gene_name']; start = rows['start']; stop =rows['stop']; chr=str(rows['chr']); pid=rows['pid']; gene_id=rows['gene_id']; gene_biotype=rows['gene_biotype']; transcript_id=rows['transcript_id']; strand=rows['strand']
            if not gene_name == rows['gene_name']:
                if gene_name is not None:
                    annotationDB.write(chr+'\t'+str(start)+'\t'+str(stop)+'\t'+str(stop-start)+'\t'+strand+'\t'+gene_name+'\t'+pid+'\t'+gene_id+'\t'+gene_biotype+'\t'+transcript_id+'\n')
                length = 0
                gene_name = rows['gene_name']; start = rows['start']; stop =rows['stop']; chr=str(rows['chr']); pid=rows['pid']; gene_id=rows['gene_id']; gene_biotype=rows['gene_biotype']; transcript_id=rows['transcript_id']; strand=rows['strand']
    annotationDB.close()
