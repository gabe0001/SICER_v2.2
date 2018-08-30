import HTSeq
import SICER_MS
import itertools
from collections import defaultdict

genome_data = SICER_MS.get_genome_data('../genomes/hg19')
iterator = HTSeq.BED_Reader('../input/sample_small_sorted_frag.bed')
window_size = 200

frag = list(itertools.islice(iterator,1))[0]


read_counts = {}

#for chrom in genome_data:
#    read_counts[chrom] = defaultdict(float)
#    window_counts_dict[chrom] = []

read_counts[frag.iv.chrom] = defaultdict(float)

# calculates number of windows that a fragment spans, then loops over each window and counts the fraction of the frag in that window. 
def pile_up_count(read_counts, frag, window_size):

    frag_start = frag.iv.start
    frag_end = frag.iv.end
    frag_size = frag_end - frag_start
    first_window_start = int(frag_start/window_size * window_size)
    first_window_end = first_window_start + window_size
    last_window_start = int(frag_end/window_size * window_size)
    last_window_end = last_window_start + window_size
    
    windows_spanned = int((last_window_end - first_window_start)/window_size)
    frag_chrom = frag.iv.chrom

    for i in range(windows_spanned):
        window_start = first_window_start + i*window_size
        window_end = first_window_end + i*window_size
        #window_iv = HTSeq.GenomicInterval(frag_chrom, window_start, window_end)
        
        if i==0:
            count_frac = float(window_end - frag_start)/float(frag_size)
            read_counts[frag_chrom][window_start] += count_frac
        elif i==(windows_spanned - 1):
            count_frac = float(frag_end - last_window_start)/float(frag_size)
            read_counts[frag_chrom][window_start] += count_frac
        #elif (i>0) and (i<windows_spanned):
        else:
            count_frac = float(window_size)/float(frag_size)
            read_counts[frag_chrom][window_start] += count_frac
