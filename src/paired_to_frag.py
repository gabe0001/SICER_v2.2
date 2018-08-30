# paired_to_frag.py
# Author: Nicholas Gabriel
# August 2018

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator
import time
import HTSeq
import bisect
import scipy
import SICER_MS


def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--bam_file", action="store", type="string", dest="bamfile", help="paired-end bam file to be condensed into a bed file of fragments")
    parser.add_option("-g", "--genome", action="store", type="string", dest="genome_data", metavar="<file>",
                      help="name of reference genome (mm9 for mouse)")
    print argv
    (opt, args) = parser.parse_args(argv)
    
    genome_file = "genomes/" + opt.genome_data

    # read genome data from file containing genome data
    # store genome data in the dictionary genome
    genome = SICER_MS.get_genome_data(genome_file)

    #Sort paired-end BAM file by name so paired reads are adjacent
    bam_sorted = opt.bamfile[:-4]+"_sorted.bam"
    os.system('samtools sort -O BAM -n %s > %s' % (opt.bamfile, bam_sorted))
    
    bam_reader = HTSeq.BAM_Reader(bam_sorted)
    
    outfile = open(opt.bamfile[:-4]+"_frag.bed", 'w')
    
    pe_iterator = HTSeq.pair_SAM_alignments(bam_reader)

    for mate1, mate2 in pe_iterator:
        # try statement for mate = None, which yields an error when mate.proper_pair is evaluated 
        try:
            if not mate1.proper_pair:
                if mate2 == None:
                    #singletons += 1
                    continue
                else:
                    #improper_pairs += 1
                    continue
            if not mate2.proper_pair:
                if mate1 == None:
                    #singletons += 1
                    continue
                else:
                    #improper_pairs += 1
                    continue
            if (mate1.iv.chrom not in genome):
                continue 
            
            if mate1.iv.start < 0:
                #improper_pairs += 1
                continue
            if mate1.iv.end > genome[mate1.iv.chrom]:
                #improper_pairs += 1
                continue
            if mate2.iv.start < 0:
                #improper_pairs += 1
                continue
            if mate2.iv.end > genome[mate2.iv.chrom]:
                #improper_pairs += 1
                continue

            elif (mate1.iv.chrom == mate2.iv.chrom) :
               
                frag_start = min([mate1.iv.start, mate2.iv.start])
                frag_end = max([mate1.iv.end, mate2.iv.end])
                
                frag_iv = HTSeq.GenomicInterval(mate1.iv.chrom, frag_start, frag_end, mate1.iv.strand)                
                                
                SICER_MS.write_to_outfile(outfile, frag_iv, mate1)
            
        except(AttributeError):
            # singleton!
            continue

   
if __name__ == "__main__":
        main(sys.argv)
