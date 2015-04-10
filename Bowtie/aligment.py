__author__ = 'peeyush'

import subprocess as sp
import os, sys
from os import listdir


# setup our program variables
program = 'bowtie2'
refgenomename = 'GRCh37'

refgenomename = "/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"


# the data directory, and read the filenames
datadir = '/home/peeyush/Desktop/Bowtie_python_test/data/samples/'


# make a new directory for our output data
outdir = '/home/peeyush/Desktop/Bowtie_python_test/data/sample_maps/'
name = 'something'


try:
    os.mkdir(outdir)
except Exception as error:
    print error
    print "WARNING:  OUTPUT DIRECTORY ALREADY EXISTS!  Data may be lost!"


# Join multiple fastq into one
newfilename = 'cat '
for i in listdir(datadir):
    if i.endswith(".gz"):
        print i
        newString = newfilename+" "+i
    newfilename = newfilename+" > "+name+".fastq.gz"

#bowtie2 -t -k 3 -p 6 -x /ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome \
#                   -U <(zcat /ps/imt/e/new_bam/H3K4me3_Non_RA/H3K4me3_Non_RA.gz) -S /ps/imt/e/new_bam/H3K4me3_Non_RA/H3K4me3_Non_RA.sam

uncompress_input = '<(zcat '+newString+')'


# use a loop to run multiple instances of the bowtie mapper
for readfn in newfilename:
    print readfn
    # make the outfile name from the readfile name, add the extension .map
    readfn = datadir + readfn
    outfn = outdir + readfn.split('/')[-1].split('.')[0] + '.sam'
    print readfn
    print outfn
    proc = sp.Popen( [program, refgenomename, readfn, outfn] )