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
        import random
        deduppath = os.path.join(self.resultdir, 'alignedLane', self.name + '_dedup', self.name + '_dedup')
        commons.ensure_path(deduppath)
        self.deduppath = os.path.join(deduppath + '_' + self.genome.name)
        bamfile = pysam.Samfile(self.bampath, "rb")
        genome = self.genome
        last_forward_position = -1
        last_reverse_position = -1
        forward_reads = set()
        reverse_reads = set()
        last_chr = -1
        count = 0
        dup = 0
        rdup = 0
        dup_dict = {}
        out_sam = pysam.Samfile('/media/peeyush/F87E9CDC7E9C94CA/M1_M2_Data/test.bam', 'wb',
                                reference_names=genome.genome.references, reference_lengths=genome.refrence_length())
        for read in bamfile.fetch():
            #print read
            #print 'chr', bamfile.references[read.tid]

            if not read.is_reverse:
                if read.pos == last_forward_position:
                    dup += 1
                    forward_reads.add(read)

                else:
                    if dup > 7:
                        dup_dict[last_forward_position]=dup
                        forward_reads = random.sample(forward_reads, 2)
                        for rd in forward_reads:
                            out_sam.write(rd)
                    forward_reads = set()
                    last_forward_position = read.pos
                    dup = 0
            else:
                readpos = read.pos + read.qlen
                if readpos == last_reverse_position:
                    reverse_reads.add(read)
                    rdup += 1

                else:
                    if rdup > 7:
                        dup_dict['r'+str(last_reverse_position)]=rdup
                        reverse_reads = random.sample(reverse_reads, 2)
                        for rr in reverse_reads:
                            out_sam.write(rr)
                    else:
                        forward_reads = set()
                        forward_reads.add(read)
                        last_reverse_position = read.pos

                    last_reverse_position = readpos
                    rdup = 0
            #print 'dup', read.opt('XC')
            last_chr = read.tid
            read_chr = bamfile.references[read.tid]
            type(read_chr)
            #filtered_on_this_chromosome = [(genome.get_chromosome_length(read_chr), genome.get_chromosome_length(read_chr))]
            #last_reverse_position = -1
            #last_forward_position = -1
            count += 1
            if count > 5000: break
        out_sam.close()
        pysam.index('/media/peeyush/F87E9CDC7E9C94CA/M1_M2_Data/test.bam')
        print dup


    def callPeaks(self, sample, controlsample, name, peakcaller):
        return