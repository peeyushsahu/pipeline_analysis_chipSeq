__author__ = 'peeyush'
from pandas import read_csv
import subprocess as sp

base_path = "/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/"
path_to_seq = "seq4motif/"
path_to_genome = "/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/"
path_to_program = "/home/peeyush/peeyush/meme/bin/"
program = "meme-chip"
motif_db_base = "/home/peeyush/peeyush/meme/db/motif_databases/"

def seq4motif(peak_data):
    '''
    This function will take dataframe with chr, start and stop position.
    :param peak_data:
    :return: Create a file of sequence for the given dataframe
    '''
    output_dir = []
    for k,df in peak_data.iteritems():
        if len(df) > 0:
            output_dir.append(k.translate(None, ' '))
            file = open(base_path+path_to_seq+k.translate(None, ' ')+".txt", "w")
            count = 1
            gc_percent = []
            CpG_ratio = []
            for v, row in df.iterrows():
                summit = int(float(row['summit'])) + int(row['start'])
                #tss = row['Next Transcript tss distance']
                start = summit - 250
                stop = summit + 250
                if 'GL' in str(row['chr']):
                    CpG_ratio.append(0)
                    gc_percent.append(0)
                    continue
                else:
                    if row['chr'] == 'X' or row['chr'] == 'Y' or len(str(row['chr'])) > 5:
                        seq = get_sequence(row['chr'], start, stop)
                        if "NNNNNN" in seq or len(seq) == 0:
                            #print row
                            print "chr", row['chr'], "Start", start, "end", stop
                    else:
                        seq = get_sequence(int(float((row['chr']))), start, stop)
                        if "NNNNNN" in seq or len(seq) == 0 :
                            print "chr", int(float((row['chr']))), "Start", start, "end", stop
                    CpG = CpG_value(seq)
                    CpG_ratio.append(CpG[0])
                    gc_percent.append(CpG[1])
                    file.write(">seq:"+str(count)+":"+row['Next transcript gene name']+"\n")
                    file.write(seq+"\n")
                count += 1
            file.close()
            print 'df',len(df)
            print 'CpG',len(CpG_ratio)
            df['CpG_ratio'] = CpG_ratio
            df['CG_percent'] = gc_percent
            df.to_csv(base_path+'/CpG/'+k+'.csv', sep=",", encoding='utf-8', ignore_index=True)
    return output_dir

def get_sequence(chr, start, end):
    '''
    This load genome for sequence fetch
    :param chr:
    :param start:
    :param end:
    :return: Sequence
    '''
    import pysam
    genome = pysam.Fastafile(path_to_genome+'genome.fa')
    sequence = genome.fetch(chr, start, end)
    return sequence

def motif_analysis(db, nmotif, output_dir):
    for seqFile in output_dir:
        output = base_path+path_to_seq+"motif_results/"+seqFile
        #database = motif_db_base+db[0]+' -db '+motif_db_base+db[1]
        input_file = base_path+path_to_seq+seqFile+'.txt'
        #print path_to_program+program, output, database, nmotifs, input_file
        proc = sp.Popen([path_to_program+program, '-o', output, '-db', motif_db_base+db[0], '-db', motif_db_base+db[1], '-meme-nmotifs', str(nmotif), input_file])
        proc.wait()


# CG% should be more than 50%
# Obs/Exp: CpG_ratio should be more then 0.6
# Obs/Exp CpG = Number of CpG / (Number of C * Number of G) * N   where N = length of sequence.
#seq = #'GCCTCCTCCGAACGCGGCCGCCTCCTCCTCCGAACGTGGCCTCCTCCGAACGCGGCCGCCTCCTCCTCCGAACGCGGCCGCCTCCTCCTCCGAACGTGGCCTCCTCCG
#AACGTGGCCGCCTCCTCCTCCGAACGTGGCCTCCTCCGAACGCGGCCGCCGCCTCCTCCGAACGCGGCCTCCTCCTCCTCCGAACGCGGCCGCCTCCTCCTCCGAACGTGGCCGCCT
#CCGAACGTGGCCGCCGCCTCCTCCGAACGTGGCCGCTTCCGCAGCGCCCGGCGCAGGCCGCACTCCGCCACCAGGGGGCGCCACAGCTCCTCGCGCCGCCGCCTCCCGCAAACACAA
#AGAGCCGCGCGGCCACGACGGCCGCGTGCCCGGAGCGCCGGGGTCTTTCCTGGGCTCCAAAGTCAAGAGCTCACGTTCCGGGAGGATCTGTCCGCGGAAATTCGGTTCTGAGCGTCG
#CCGGACTCCGCCGCGGGGAGGCGGGTGAGGGGAGGGGGCCG'

def CpG_value(seq):
    '''
    This function will calculate CpG enrichment for the sequence
    :param seq:
    :return:
    '''
    seq = seq.upper()
    return float(seq.count('CG')) / (seq.count('C') * seq.count('G') ) * len(seq), float(100)/len(seq) * (seq.count('C') + seq.count('G'))

#peaks_list = ['Unique_PRMT6_2_RA_seq6 vs IgG_RA_seq6 filtered', 'Unique_PRMT6_2_seq6 vs IgG_seq6 filtered']
#peak_data = {}
#for a in peaks_list:
#    df = read_csv(
#        '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/overlap/non_overlapping/' + a + '.csv',
#        header=0, sep=',')
#    peak_data[a] = df
#seq4motif(peak_data)
