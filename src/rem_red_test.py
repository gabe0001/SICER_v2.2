import SICER_MS
import HTSeq

test_type = 'bam'

if test_type == 'bed':
    infile = '../input/sample_small_sorted.bed'
    outfile = '../input/rem_red_test.bed'
    genome_data = SICER_MS.get_genome_data('../genomes/hg19')
    cutoff = 1
    total, retained = SICER_MS.remove_redundant_reads_bed(infile, outfile, cutoff, genome_data)
    
if test_type == 'bam':
    infile = '../input/sample_small_sorted.bam'
    outfile = '../input/rem_red_test.bam'
    genome_data = SICER_MS.get_genome_data('../genomes/hg19')
    cutoff = 1
    total, retained = SICER_MS.remove_redundant_reads_bam(infile, outfile, cutoff, genome_data)


