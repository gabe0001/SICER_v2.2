# Author: Ujwal Boddeti
# Version 1.0  11/10/2017

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator
import time
import HTSeq
import bisect
import scipy

def write_to_outfile(outfile, GI):
    outfile.write(str(GI.chrom)+"\t"+str(GI.start)+"\t"+str(GI.end)+"\t"+str(GI.strand)+"\n")

def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--bam_file", action="store", type="string", dest="bamfile", help="paired-end bam file to be condensed into a single-ended bam file")
    print argv
    (opt, args) = parser.parse_args(argv)

    #Sortng inputted paired-end BAM file by name, so paired reads are adjacent
    bam_reader = HTSeq.BAM_Reader(opt.bamfile)
    os.system('samtools sort -O BAM -n %s > %s' % (opt.bamfile, opt.bamfile[:-4]+"_sorted.bam"))

    #BAM file is converted to BED format
    os.system('bamToBed -i %s > %s' % (opt.bamfile[:-4]+"_sorted.bam", opt.bamfile[:-4]+"_sorted.bed"))

    #BED iterator created to traverse BED file
    bed_iterator = HTSeq.BED_Reader(opt.bamfile[:-4]+"_sorted.bed")

    outfile = open(opt.bamfile[:-4]+"_sorted_condensed.bed", 'w')

    #algorithmic logic -> combine paired reads into a single read
    pair1 = None
    pair2 = None
    oddRead = None
    new_start = 0
    new_end = 0
    pos = 0
    new_strand = ''
    singleRead_iv = None
    line_count = 0
    error_count = 0
    finalcount = 0

    for read in bed_iterator:
        line_count =+ line_count + 1
        if line_count%2 != 0:
            oddRead = read.iv

        elif line_count%2 == 0:
            pair1 = oddRead
            pair2 = read.iv
            #print "Read 1 chr: " + str(pair1.chrom) + " and Read 2 chr: " + str(pair2.chrom)

            if str(pair1.chrom) == str(pair2.chrom):
                #determines start and end of single read
                new_start = min([pair1.start, pair2.start])
                new_end = max([pair1.end, pair2.end])

                #position calculation
                pos = (new_start + new_end) / 2
                if pos%1 != 0:
                    pos=pos-0.5

                #decides proper strand for new single read (strand is that of leftmost read in pair)
                if pair1.start < pair2.start:
                    new_strand = pair1.strand
                else:
                    new_strand = pair2.strand

                #creates new single read with length 1 bp
                singleRead_iv = HTSeq.GenomicInterval(pair1.chrom, pos, pos+1, new_strand)

                #writes read to output BED file
                write_to_outfile(outfile, singleRead_iv)
                finalcount = finalcount + 1
            else:
                #print "Error: paired reads not on same chromosome."
                error_count = error_count + 1



    #print "pairedRead to singleRead conversion is done running.\nThere were " + str(error_count) + " pairs on different chromosomes"
    #print "Reads written to outfile: " + str(finalcount)

if __name__ == "__main__":
        main(sys.argv)
