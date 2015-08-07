from alignment import aligner

__author__ = 'peeyush'
import subprocess as sp
import os
import timeit
import commons

class Lane():


    def __init__(self, samplename, path):
        self.name = samplename
        self.path = path
        self.genome = None
        self.fqpath = None
        self.resultdir = None
        self.sampath = None


    def join_multiple_fq(self):
        self.resultdir = commons.create_odir()
        fqdir = os.path.join(self.resultdir, 'cache', self.name)
        commons.create_sample_dir(fqdir)
        path = self.path
        outname = 'temp_' + self.name + '.fastq.gz'
        self.fqoutpath = os.path.join(fqdir, outname)
        newfilename = 'cat '
        for i in os.listdir(path):
            if i.endswith("fastq.gz"):
                newfilename = newfilename+" "+path+i
        newfilename = newfilename+" > "+self.fqoutpath
        print newfilename
        proc = sp.Popen([newfilename], shell=True)
        proc.wait()
        #os.remove(self.fqoutpath)
        #print outname


    def do_alignment(self, genome):
        aligner.bowtie2(self, genome)
        return


    def sam2bam(self):
        import pysam
        return

class AlignedLane():


    def __init__(self, Lane):
        self.name = Lane.name
        self.path = Lane.path
        self.sampath = None
        self.reultpath = Lane.resultdir
