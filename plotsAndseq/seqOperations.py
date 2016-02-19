__author__ = 'peeyush'
from pandas import read_csv
import subprocess as sp

base_path = "/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/"
path_to_seq = "seq4motif/"
path_to_genome = "/ps/imt/f/Genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/"
path_to_program = "/home/sahu/meme/bin/"
program = "meme-chip"
motif_db_base = "/home/sahu/Documents/motif_databases/"

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
                start = summit - 500
                stop = summit + 500
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
            df = df[df['CpG_ratio'] >= 0.6]
            df = df[df['CG_percent'] >= 50]
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
    '''
    This method calls meme_chip
    :param db:
    :param nmotif:
    :param output_dir:
    :return:
    '''
    for seqFile in output_dir:
        output = base_path+path_to_seq+"motif_results/"+seqFile
        input_file = base_path+path_to_seq+seqFile+'.txt'
        meme_chip(output, input_file, db, nmotif)



def meme_chip(output, input_file, db, nmotif):
    '''
    Initialize meme-chip on given seq file.
    :param output: path to store output
    :param input_file: path to the seqfile
    :param db: databases used for motif enrichment
    :param nmotif: No. of motifs to find minimum
    :return:
    '''
    proc = sp.Popen([path_to_program+program, '-noecho', '-o', output, '-db', motif_db_base+db[0], '-db', motif_db_base+db[1], '-meme-nmotifs', str(nmotif), input_file])
    proc.wait()

def centrimo(output, input_file, neg_control, db):
    '''
    Runs centrimo analysis to obtain differential enrichment of motifs in input_file compared to neg_control
    :param output: path to store output
    :param input_file: path to the seqfile
    :param neg_control: path to the seqfile
    :param db: databases used for motif enrichment
    :return:
    '''
    proc = sp.Popen([path_to_program+'centrimo', '-score', 5.0, '-ethresh', 10.0, '-neg', neg_control,'-verbosity', 1, '-o', output, '-db', motif_db_base+db[0], '-db', motif_db_base+db[1], input_file])
    proc.wait()

def CpG_value(seq):
    '''
    This function will calculate CpG enrichment for the sequence.
    # CG% should be more than 50%
    # Obs/Exp: CpG_ratio should be more then 0.6
    # Obs/Exp CpG = Number of CpG / (Number of C * Number of G) * N   where N = length of sequence.
    :param seq:
    :return:
    '''
    seq = seq.upper()
    return float(seq.count('CG')) / (seq.count('C') * seq.count('G') ) * len(seq), float(100)/len(seq) * (seq.count('C') + seq.count('G'))


def density_based_motif_comparision(dataframe, columnname):
    '''
    This method will take a dataframe and sort it according to columnname and divide into four parts.
    Sequences for the peaks for every part extracted and saved. 0: low peak strength; 3: high peak strength.
    :param dataframe:
    :param columnname:
    :return:
    '''
    df = dataframe
    print "Density based motif analysis: Sorting by column "+columnname
    df = df.sort(columnname, ascending=True)
    size = df.shape[0]
    count = 0
    dfstart = 0
    dfstop = size/4
    for i in range(0,4):
        file = open(base_path+'density_based_motif/'+columnname.translate(None, ' ')+str(i)+".txt", "w")
        for v, row in df[dfstart:dfstop].iterrows():
            summit = int(float(row['summit'])) + int(row['start'])
            #tss = row['Next Transcript tss distance']
            start = summit - 250
            stop = summit + 250
            if 'GL' in str(row['chr']):
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
                file.write(">seq:"+str(count)+":"+row['Next transcript gene name']+"\n")
                file.write(seq+"\n")
            count += 1
        dfstart = dfstop
        dfstop = dfstart+size/4
        file.close()
        print "Motif search for:"+columnname+str(i)
        meme_chip(base_path+'density_based_motif/'+columnname.translate(None, ' ')+str(i), base_path+'density_based_motif/'+columnname.translate(None, ' ')+str(i)+".txt",
                  ["JASPAR_CORE_2014_vertebrates.meme", "uniprobe_mouse.meme"], 10)
        print dfstart, dfstop