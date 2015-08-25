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
        coverage_before = []
        coverage_after = []
        dup_dict = {}
        last_forward_position = -1
        last_reverse_position = -1
        forward_reads = set()
        reverse_reads = set()
        last_chr = -1
        last_read = 0
        count = 0
        out_sam = pysam.Samfile(self.deduppath+'.temp', 'wb',
                                reference_names=genome.genome.references, reference_lengths=genome.refrence_length())
        last_chr = -1
        for read in bamfile.fetch():
            #print read
            #print 'chr', bamfile.references[read.tid]
            count += 1
            if not last_chr == read.tid:
                last_chr = read.tid
                read_chr = bamfile.references[read.tid]
                print 'chr Start', bamfile.references[read.tid], 'read_pos', read.pos, 'Read count', count

            if not read.is_reverse:
                if read.pos == last_forward_position:
                    forward_reads.add(read)

                else:
                    dup_dict['+'+str(len(forward_reads))] = dup_dict.get('+'+str(len(forward_reads)), 0) + 1
                    if len(forward_reads) < 7:
                        if len(forward_reads) >= 2:
                            forward_reads = random.sample(forward_reads, 2)
                        for rd in forward_reads:
                            out_sam.write(rd)
                    if len(forward_reads) > 7:
                        forward_reads = random.sample(forward_reads, 1)
                        out_sam.write(forward_reads.pop())
                    forward_reads = set()
                    forward_reads.add(read)
                    last_forward_position = read.pos
            else:
                readpos = read.pos + read.qlen
                if readpos == last_reverse_position:
                    reverse_reads.add(read)

                else:
                    dup_dict['-'+str(len(forward_reads))] = dup_dict.get('-'+str(len(forward_reads)), 0) + 1
                    if len(reverse_reads) < 7:
                        if len(reverse_reads) > 2:
                            reverse_reads = random.sample(reverse_reads, 2)
                        for rr in reverse_reads:
                            out_sam.write(rr)
                    if len(reverse_reads) > 7:
                        reverse_reads = random.sample(reverse_reads, 1)
                        out_sam.write(reverse_reads.pop())
                    reverse_reads = set()
                    reverse_reads.add(read)
                    last_reverse_position = readpos

        # Last push for reads
        if len(forward_reads) < 7:
            if len(forward_reads) >= 2:
                forward_reads = random.sample(forward_reads, 2)
            for rd in forward_reads:
                out_sam.write(rd)
        if len(forward_reads) > 7:
            forward_reads = random.sample(forward_reads, 1)
            out_sam.write(forward_reads.pop())
        if len(reverse_reads) < 7:
            if len(reverse_reads) > 2:
                reverse_reads = random.sample(reverse_reads, 2)
            for rr in reverse_reads:
                out_sam.write(rr)
        if len(reverse_reads) > 7:
            reverse_reads = random.sample(reverse_reads, 1)
            out_sam.write(reverse_reads.pop())

            #print 'dup', read.opt('XC')
            #filtered_on_this_chromosome = [(genome.get_chromosome_length(read_chr), genome.get_chromosome_length(read_chr))]
            #if count > 15000: break
        out_sam.close()
        pysam.sort(self.deduppath+'.temp', self.deduppath)
        pysam.index(self.deduppath+'.bam')
        #print dup


    def callPeaks(self, sample, controlsample, name, peakcaller):
        return