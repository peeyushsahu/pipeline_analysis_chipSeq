from alignment import aligner

__author__ = 'peeyush'
import subprocess as sp
import os
import timeit
import commons
import pysam
import random


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
                newfilename = os.path.join(newfilename+" "+path,i)
        newfilename = newfilename+" > "+self.fqoutpath
        print newfilename
        proc = sp.Popen([newfilename], shell=True)
        proc.wait()
        #os.remove(self.fqoutpath)
        print outname


    def do_alignment(self, genome, method):
        self.genome = genome
        if method == "Bowtie2":
            aligner.bowtie2_aligner(self, genome)
            self.sam2bam()
            self.bam_sort()
        if method == "Tophat2":
            aligner.tophat2_aligner(self, genome)
            self.bam_sort()
        return


    def sam2bam(self):
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


    def bam_sort(self):
        self.sortbampath = os.path.join(self.resultdir, 'alignedLane', self.name, self.name + '_' + self.genome.name)
        print self.sortbampath
        pysam.sort(self.bampath, self.sortbampath)
        self.bampath = self.sortbampath+'.bam'
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

    def do_dedup(self, maximum_stacks=None, maximum_stacks_allowed=2):
        """
        To remove PCR duplicates from bam files.
        """
        deduppath = os.path.join(self.resultdir, 'alignedLane', self.name + '_dedup', self.name + '_dedup')
        commons.ensure_path(deduppath)
        self.deduppath = os.path.join(deduppath + '_' + self.genome.name)
        bamfile = pysam.Samfile(self.bampath, "rb")
        genome = self.genome
        dup_dict = {}
        last_forward_position = -1
        last_reverse_position = -1
        forward_reads = set()
        reverse_reads = set()
        out_sam = pysam.Samfile(self.deduppath+'.temp', 'wb',
                                reference_names=genome.genome.references, reference_lengths=genome.refrence_length())
        for read in bamfile.fetch():
            if not read.is_reverse:
                if read.pos == last_forward_position:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in xrange(0,repeat_count):
                        forward_reads.add(read)
                else:
                    dup_dict['+'+str(len(forward_reads))] = dup_dict.get('+'+str(len(forward_reads)), 0) + 1
                    if maximum_stacks is None or len(forward_reads) < maximum_stacks:
                        if len(forward_reads) >= maximum_stacks_allowed:
                            forward_reads = random.sample(forward_reads, 2)
                        for rd in forward_reads:
                            out_sam.write(rd)
                    else:
                        forward_reads = random.sample(forward_reads, 1)
                        out_sam.write(forward_reads.pop())
                    forward_reads = set()
                    forward_reads.add(read)
                    last_forward_position = read.pos
            else:
                readpos = read.pos + read.qlen
                if readpos == last_reverse_position:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in xrange(0, repeat_count):
                        reverse_reads.add(read)
                else:
                    dup_dict['-'+str(len(reverse_reads))] = dup_dict.get('-'+str(len(reverse_reads)), 0) + 1
                    if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
                        if len(reverse_reads) >= maximum_stacks_allowed:
                            reverse_reads = random.sample(reverse_reads, 2)
                        for rr in reverse_reads:
                            out_sam.write(rr)
                    else:
                        reverse_reads = random.sample(reverse_reads, 1)
                        out_sam.write(reverse_reads.pop())
                    reverse_reads = set()
                    reverse_reads.add(read)
                    last_reverse_position = readpos
        # Last push for reads
        if maximum_stacks is None or len(forward_reads) < maximum_stacks:
            if len(forward_reads) >= maximum_stacks_allowed:
                forward_reads = random.sample(forward_reads, 2)
            for rd in forward_reads:
                out_sam.write(rd)
        else:
            forward_reads = random.sample(forward_reads, 1)
            out_sam.write(forward_reads.pop())
        if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
            if len(reverse_reads) >= maximum_stacks_allowed:
                reverse_reads = random.sample(reverse_reads, 2)
            for rr in reverse_reads:
                out_sam.write(rr)
        else:
            reverse_reads = random.sample(reverse_reads, 1)
            out_sam.write(reverse_reads.pop())
        out_sam.close()
        pysam.sort(self.deduppath+'.temp', self.deduppath)
        pysam.index(self.deduppath+'.bam')

    def do_dedup_stringent(self, maximum_stacks=None, maximum_stacks_allowed=2):
        """
        To remove PCR duplicates from bam files.
        """
        deduppath = os.path.join(self.resultdir, 'alignedLane', self.name + '_dedup', self.name + '_dedup')
        commons.ensure_path(deduppath)
        self.deduppath = os.path.join(deduppath + '_' + self.genome.name)
        bamfile = pysam.Samfile(self.bampath, "rb")
        genome = self.genome
        dup_dict = {}
        last_forward_start = -1
        last_forward_end = -1
        last_reverse_start = -1
        last_reverse_end = -1
        forward_reads = set()
        reverse_reads = set()
        out_sam = pysam.Samfile(self.deduppath+'.temp', 'wb',
                                reference_names=genome.genome.references, reference_lengths=genome.refrence_length())
        for read in bamfile.fetch():
            if not read.is_reverse:
                if read.pos == last_forward_start and read.pos+read.qlen == last_forward_end:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in xrange(0,repeat_count):
                        forward_reads.add(read)
                else:
                    dup_dict['+'+str(len(forward_reads))] = dup_dict.get('+'+str(len(forward_reads)), 0) + 1
                    if maximum_stacks is None or len(forward_reads) < maximum_stacks:
                        if len(forward_reads) >= maximum_stacks_allowed:
                            forward_reads = random.sample(forward_reads, 2)
                        for rd in forward_reads:
                            out_sam.write(rd)
                    else:
                        forward_reads = random.sample(forward_reads, 1)
                        out_sam.write(forward_reads.pop())
                    forward_reads = set()
                    forward_reads.add(read)
                    last_forward_start = read.pos
                    last_forward_end = read.pos+read.qlen
            else:
                readpos = read.pos + read.qlen
                if readpos == last_reverse_start and read.pos == last_reverse_end:
                    try:
                        repeat_count = read.opt('XC')
                    except KeyError: #no XC
                        repeat_count = 1
                    for ii in xrange(0, repeat_count):
                        reverse_reads.add(read)
                else:
                    dup_dict['-'+str(len(reverse_reads))] = dup_dict.get('-'+str(len(reverse_reads)), 0) + 1
                    if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
                        if len(reverse_reads) >= maximum_stacks_allowed:
                            reverse_reads = random.sample(reverse_reads, 2)
                        for rr in reverse_reads:
                            out_sam.write(rr)
                    else:
                        reverse_reads = random.sample(reverse_reads, 1)
                        out_sam.write(reverse_reads.pop())
                    reverse_reads = set()
                    reverse_reads.add(read)
                    last_reverse_end = read.pos
                    last_reverse_start = readpos
        # Last push for reads
        if maximum_stacks is None or len(forward_reads) < maximum_stacks:
            if len(forward_reads) >= maximum_stacks_allowed:
                forward_reads = random.sample(forward_reads, 2)
            for rd in forward_reads:
                out_sam.write(rd)
        else:
            forward_reads = random.sample(forward_reads, 1)
            out_sam.write(forward_reads.pop())
        if maximum_stacks is None or len(reverse_reads) < maximum_stacks:
            if len(reverse_reads) >= maximum_stacks_allowed:
                reverse_reads = random.sample(reverse_reads, 2)
            for rr in reverse_reads:
                out_sam.write(rr)
        else:
            reverse_reads = random.sample(reverse_reads, 1)
            out_sam.write(reverse_reads.pop())
        out_sam.close()
        pysam.sort(self.deduppath+'.temp', self.deduppath)
        pysam.index(self.deduppath+'.bam')

    def callPeaks(self, peakscaller, sample, controlsample, name, outdir, broad_cutoff, broadpeaks=False):
        return