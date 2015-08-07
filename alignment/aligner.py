__author__ = 'peeyush'

import subprocess as sp
import os, sys
import commons


class human_GRCh37():
    def __init__(self):
        self.name = 'GRCh37'
        self.refgenomename = "/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"


def bowtie2(lane, genome):
    # setup our program variables
    print 'Bowtie alignment'
    program = '/home/peeyush/Documents/bowtie2-2.2.4/bowtie2'
    # make the outfile name from the readfile name, add the extension .map
    thread = '-p 6'
    unaligned = '--un ' + os.path.join(lane.resultdir, 'cache', lane.name, lane.name+'_un_aligned.fastq')
    readfn = lane.fqoutpath
    outpath = os.path.join(lane.resultdir, 'alignedLane', lane.name)
    commons.create_sample_dir(outpath)
    outfn = os.path.join(outpath, lane.name + '_' + genome.name + '.sam')
    lane.sampath = outfn
    filename = os.path.join(outpath, lane.name + '_' + genome.name + '_bowtie_stats.txt')
    stats = open(filename, 'wb')
    cmd = ' '.join([program, thread, unaligned, '-x', genome.refgenomename, '-U', readfn, '-S', outfn])
    print cmd
    proc = sp.Popen(cmd, stdout=stats, shell=True)
    proc.wait()
    stats.flush()