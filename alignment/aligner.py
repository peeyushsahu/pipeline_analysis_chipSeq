__author__ = 'peeyush'

import subprocess as sp
import os, sys
import commons


class human_GRCh37():
    def __init__(self):
        self.name = 'GRCh37'
        self.refgenomename = "/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"


class samtool():
    def __init__(self):
        self.samtools = '/home/peeyush/Documents/samtools-1.1/samtools'

def sample_dir(lane):
    a = ['cache', 'peaks', 'alignedLane']
    for i in a:
        commons.ensure_path(os.path.join(lane.resultdir, i, lane.name))


def bowtie2(lane, genome):
    # setup our program variables
    sample_dir(lane)
    print 'Bowtie alignment'
    program = '/home/peeyush/Documents/bowtie2-2.2.4/bowtie2'
    # make the outfile name from the readfile name, add the extension .map
    thread = '-p 6'
    unaligned = '--un ' + os.path.join(lane.resultdir, 'cache', lane.name, lane.name+'_un_aligned.fastq')
    readfn = lane.fqoutpath
    outpath = os.path.join(lane.resultdir, 'cache', lane.name)
    lane.sampath = os.path.join(outpath, lane.name + '_' + genome.name + '.sam')
    lane.bampath = os.path.join(outpath, lane.name + '_' + genome.name + '.bam')
    lane.temp_files.append(lane.sampath, lane.bampath)
    statfile = os.path.join(lane.resultdir, 'alignedLane', lane.name, lane.name + '_' + genome.name + '_bowtie_stats.txt')
    cmd = ' '.join([program, thread, unaligned, '-x', genome.refgenomename, '-U', readfn, '-S', lane.sampath])
    print cmd
    try:
     proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
     stdout, stdrr = proc.communicate()
     with open(statfile, 'wb') as op:
         op.write(stdrr+"\n") ##Writing bowtie stats
     op.close()
     proc.wait()
    except:
        raise IOError ('Subprocess Bowtie2 exited with error:',proc)