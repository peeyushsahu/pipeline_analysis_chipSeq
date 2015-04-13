__author__ = 'peeyush'

def parse_gtf(path):
    from pandas import read_csv
    colnames = ['chr', 'transcript_annotation', 'region', 'start', 'stop', 'score', 'strand', 'frame', 'further_information']
    gtf_file = read_csv(path, sep='\t', colnames=0)
    gtf_file.columns = colnames
    return gtf_file

def exon_intron_vector(gtf_file):
    import numpy as np
    genome_vector = np.zeros(1565134, dtype=np.uint32)
    region = ['exon', 'transcript']
    for i in region:
        if i == 'exon':
            for k,v in gtf_file[gtf_file['region'] == i].iterrows():
                genome_vector[v['start']: v['stop']] = 2

        if i == 'transcript':
            for k,v in gtf_file[gtf_file['region'] == i].iterrows():
                genome_vector[v['start']] = 1
    return genome_vector

# unique_elements = set(genome_vector[start:stop])
# unique_count = np.bincount(genome_vector[start:stop])
# zip_unique = zip(unique_elements, unique_count)

def annotate_IntronExon_junction(dataFrame, genome_vector):
