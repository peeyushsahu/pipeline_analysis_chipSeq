__author__ = 'peeyush'



def parse_gtf(path):
    from pandas import read_csv
    #gtf_file = read_csv('/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf', sep='\t', header=0)
    colnames = ['chr', 'transcript_annotation', 'region', 'start', 'stop', 'score', 'strand', 'frame', 'further_information']
    gtf_file = read_csv(path, sep='\t')
    gtf_file.columns = colnames
    return gtf_file

def vectors_4_chromosome():
    import numpy as np
    chr = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    position = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,
            133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566, 155270560, 59373566]
    chr_len = dict(zip(chr, position))
    chr_vector = {}
    for k, v in chr_len.iteritems():
        chr_vector[k] = np.zeros(v, dtype=np.uint32)
    return chr_vector

def exon_intron_vector(path):
    import re
    print "reading GTF file"
    gtf_file = parse_gtf(path)
    print "creating chromosome vector"
    chr_vector = vectors_4_chromosome()
    print "Annotating chr_vector"
    region = ['exon', 'transcript']
    for i in region:
        if i == 'exon':
            print 'exon'
            gene_name = ''
            index = 0
            stop = 0
            for k, v in gtf_file[gtf_file['region'] == i].iterrows():
                if len(str(v['chr'])) < 5 and v['chr'] != 'MT':
                    chr_vector[str(v['chr'])][v['start']: v['stop']] = 2
                    #print v['start'], v['stop']
                    name = re.split(';| |"', v['further_information']).index('gene_name')
                    new_gene_name = re.split(';| |"', v['further_information'])[name+2]
                    if index > 0 and new_gene_name == gene_name:
                    #print gene_name
                    #print index
                        chr_vector[str(v['chr'])][stop+1: v['start']-1] = 3
                        gene_name = new_gene_name
                        index += 1
                    stop = v['stop']
                    gene_name = new_gene_name
                    index += 1
        if i == 'transcript':
            print 'transcript'
            for k, v in gtf_file[gtf_file['region'] == i].iterrows():
                if len(str(v['chr'])) < 5 and v['chr'] != 'MT':
                    chr_vector[str(v['chr'])][v['start']] = 1
                    if str(v['chr'])=='X':
                        print k, 'Start', v['start']
    del(gtf_file)
    return chr_vector

def annotate_intronexon_junction(df, chr_vector):
    '''
    Calculate Intron-Exon juction will be called for the cookie-cut (+-250) peaks from summit
    :param df: Dataframe
    :param chr_vector: annotated chromosome vector
    :return:
    '''
    import numpy as np
    dataFrame = df
    dataFrame[['chr']] = dataFrame[['chr']].astype(str)
    dataFrame['exon_intron_junction'] = 0
    dataFrame['exon_intron_percent'] = 0
    for k, v in dataFrame.iterrows():
        if len(v['chr']) < 5:
            #print v
            #print type(v['chr']),type(v['start']),type(v['stop'])
            summit = v['start']+v['summit']
            unique_elements = set(chr_vector[v['chr']][summit-250:summit+250])
            unique_count = np.bincount(chr_vector[v['chr']][v['start']:v['stop']])
            zip_unique = dict(zip(unique_elements, unique_count))
            if unique_elements.issuperset([2,3]) and zip_unique[2]>0 and zip_unique[3]>0:
                if unique_elements.issuperset([0]) and zip_unique[0]<10 or unique_elements.issuperset([0]) == False:
                    if float(zip_unique[3])/float(zip_unique[2]) < 2 and float(zip_unique[3]) / float(zip_unique[2]) > 0.5:
                        print zip_unique
                        dataFrame.loc[k, 'exon_intron_junction'] = 1
                        dataFrame.loc[k, 'exon_intron_percent'] = float(zip_unique[3])/float(zip_unique[2])
    dataFrame.to_csv('/home/peeyush/Desktop/exon_intron_junction_PRMT6_RA.csv', sep=",", encoding='utf-8', ignore_index=True, index=False)
    return dataFrame
