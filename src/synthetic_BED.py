## synthetic_BED.py
## Nick Gabriel
## June 2018

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator
import time
import HTSeq
import bisect
import scipy
import numpy as np
import SICER_MS 

def synthetic_bed(genome_path, num_reads):
		
		outfile = open('synthetic_' + genome_path.split('/')[-1] + '_data' + '.bed', 'w')
		genome = SICER_MS.get_genome_data(genome_path)
		num_chroms = len(genome.keys())
		max_read_len = 1000
		strands = ['+','-']
		
		for line in range(num_reads):
		
			chrom_idx = np.random.randint(0, num_chroms)
			rand_chrom = genome.keys()[chrom_idx]
			rand_start = np.random.randint(0, genome[rand_chrom])
			rand_end = rand_start + np.random.randint(1, max_read_len)
			rand_strand = strands[np.random.randint(0,2)]
			
			read_out = str(rand_chrom)  + '\t' + str(rand_start) + '\t' + str(rand_end) + '\t' + str(rand_strand) + '\n'
			outfile.write(read_out)

			
