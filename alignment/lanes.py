from alignment import aligner

__author__ = 'peeyush'
import subprocess as sp
import os
import timeit
import commons
import pysam

class Lane():


    def __init__(self, samplename, path):
        self.name = samplename
        self.path = path
        self.genome = None
        self.fqpath = None
        self.resultdir = None
        self.sampath = None
        self.bampath = None
        self.temp_files = []


    def join_multiple_fq(self):
        self.resultdir = commons.create_odir()
        fqdir = os.path.join(self.resultdir, 'cache', self.name)
        commons.ensure_path(fqdir)
        path = self.path
        outname = 'temp_' + self.name + '.fastq.gz'
        self.fqoutpath = os.path.join(fqdir, outname)
        self.temp_files.append(self.fqoutpath)
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
        self.genome = genome
        aligner.bowtie2(self, genome)
        self.sam2bam()
        return


    def sam2bam(self):
        import pysam
        #samtools view -Sb alignment_rep_prmt6+.sam > alignment_rep_PRMT6+.bam
        samtool = aligner.samtool()
        samtools = samtool.samtools
        print samtools, 'view -Sb', self.sampath, '>', self.bampath
        cmd = ' '.join([samtools, 'view -Sb', self.sampath, '>', self.bampath])
        try:
            proc = sp.Popen(cmd, shell=True)
            proc.wait()
        except:
            raise IOError("Problem with samtools sam 2 bam.")
        sortpath = os.path.join(self.resultdir, 'alignedLane', self.name, self.name + '_' + self.genome.name)
        print sortpath
        pysam.sort(self.bampath, sortpath)
        self.bampath = sortpath+'.bam'
        pysam.index(self.bampath)
        self.remove_temp()

    def remove_temp(self):
        for i in self.temp_files:
            os.remove(i)


class AlignedLaneDedup():
    def __init__(self, Lane):
        self.name = Lane.name
        self.genome = Lane.genome
        self.bampath = Lane.bampath
        self.deduppath = None
        self.resultdir = Lane.resultdir
        self.peakdata = None


    def do_dedup(self):
        deduppath = os.path.join(self.resultdir, 'alignedLane', self.name + '_dedup', self.name + '_dedup')
        commons.ensure_path(deduppath)
        self.deduppath = os.path.join(deduppath + '_' + self.genome.name)
        bamfile = pysam.Samfile(self.bampath, "rb")
        genome = self.genome
        last_forward_position = -1
        last_reverse_position = -1
        last_chr = -1
        count = 0
        for read in bamfile.fetch():
            a = read
            #print read
            print 'chr', bamfile.references[read.tid]
            print 'position', read.pos
            print 'reverse', read.is_reverse
            print 'dup', read.opt('XC')
            last_chr = read.tid
            read_chr = bamfile.references[read.tid]
            type(read_chr)
            filtered_on_this_chromosome = [(genome.get_chromosome_length(read_chr), genome.get_chromosome_length(read_chr))]
            last_reverse_position = -1
            last_forward_position = -1
            count += 1
            if count > 20: break


    def callPeaks(self, sample, controlsample, name, peakcaller):
        return