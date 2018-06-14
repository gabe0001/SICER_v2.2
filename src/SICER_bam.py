# SICER_bam

# Version: 2.2
# Changes include modularization of essentials functions in SICER_MS script

# Version: 2.1
# Changes include rewriting the remove_redundant_reads_bam module which works to
# remove unmapped reads, reads belonging to chrM, and duplicate reads

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator
import time
import Background_island_probscore_statistics
import HTSeq
import bisect
import scipy
import scipy.stats
import SICER_MS

def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--file", action="store", type="string", dest="file_name", metavar="<file>",
                      help="name of BAM file")
    parser.add_option("-c", "--control", action="store", type="string", dest="control_file_name", metavar="<file>",
                      help="name of control BAM file")
    parser.add_option("-g", "--genome", action="store", type="string", dest="genome_data", metavar="<file>",
                      help="name of reference genome (mm9 for mouse)")
    parser.add_option("-r", "--redundancy", action="store", type="int", dest="redundancy",
                      metavar="<int>", help="redundancy threshold")
    parser.add_option("-w", "--window_size", action="store", type="int", dest="window_size", metavar="<int>",
                      help="size of windows used to partition genome (200 for histones, 50 for TFs")
    parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>",
                      help="fragment size determines the shift (half of fragment_size of ChIP-seq read position, in bps)")
    parser.add_option("-p", "--genome_fraction", action="store", type="float", dest="genome_fraction", metavar="<int>",
                      help="effective genome fraction: 0.8 in most cases")
    parser.add_option("-s", "--gap_size", action="store", type="int", dest="gap_size", metavar="<int>",
                      help="maximum number of base pairs between windows in the same island (usually same as window size)")
    parser.add_option("-d", "--FDR", action="store", type="string", dest="FDR", metavar="<string>",
                      help="false discovery rate controlling significance")
    parser.add_option("-i", "--input_dir", action="store", type="string", dest="input_dir", metavar="<string>",
                      help="path to input directory")
    parser.add_option("-o", "--output_dir", action="store", type="string", dest="output_dir", metavar="<string>",
                      help="path to output directory")
    parser.add_option("-a", "--SICER_dir", action="store", type="string", dest="sicer_dir", metavar="<string>",
                      help="path to directory containing SICER files")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)

    # create string names for files
    file_name = opt.file_name[:-4]
    control_file_name = opt.control_file_name[:-4]
    sorted_bam_file_name = file_name + "_sorted.bam"
    red_rem_bam_file_name = file_name + "-" + str(opt.redundancy) + "-removed.bam"
    island_bed_file_name = file_name + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + "-FDR" + str(opt.FDR) + "-island.bed"
    sorted_control_file_name = control_file_name + "_sorted.bam"
    red_rem_control_file_name = control_file_name + "-" + str(opt.redundancy) + "-removed.bam"
    # This file stores the sample summary graph in bedgraph format
    normalized_bedgraph_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-normalized.bedgraph"
    # This file stores the control summary graph in bedgraph format
    control_normalized_bedgraph_file_name = opt.output_dir + "/" + control_file_name + "-W" + str(opt.window_size) + "-normalized.bedgraph"
    # This file stores the candidate islands.
    score_island_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + ".scoreisland"
    # These files store the summary graphs.
    graph_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + ".graph"
    control_graph_file_name = opt.output_dir + "/" + control_file_name + "-W" + str(opt.window_size) + ".graph"
    # This file stores the island-filtered non-redundant raw reads
    island_filtered_file_name = opt.output_dir + "/" + file_name  + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + "-FDR" + str(opt.FDR) + "-islandfiltered.bed"
    # This file stores normalized summary graph made by the island-filtered non-redundant raw reads in bedgraph format
    islandfiltered_normalized_bedgraph_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + "-FDR" + str(opt.FDR) + "-islandfiltered-normalized.bedgraph"
    # This file stores the summary of candidate islands, including chrom start end read-count_sample read-count-control pvalue, fold change and qvalue
    islandsummary_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + "-islands-summary"
    # This file stores the summary of significant islands identified with FDR criterion.
    filtered_island_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + "-islands-summary-FDR" + str(opt.FDR)

    # convert FDR to float
    FDR = float(opt.FDR)

    genome_file = opt.sicer_dir + "/genomes/" + opt.genome_data

    # read genome data from file containing genome data
    # store genome data in the dictionary genome
    genome = SICER_MS.get_genome_data(genome_file)

    # number of islands expected in random background.
    # the E value is used for identification of candidate islands that exhibit clustering.
    e_value = 1000

    # sort bam file by chromosome then coordinate using samtools
    print "\nSorting BAM file..."
    os.system('samtools sort -O BAM %s > %s' % (opt.file_name, sorted_bam_file_name))
    # sort control bam file by chromosome then coordinate using samtools
    print "Sorting control BAM file..."
    os.system('samtools sort -O BAM %s > %s' % (opt.control_file_name, sorted_control_file_name))


    # remove redundant reads in sorted BAM file and count number of total reads and number of retained reads
    print "\nPreprocess the sorted BAM file to remove redundancy with threshold " + str(opt.redundancy) + "..."
    total, retained = SICER_MS.remove_redundant_reads_bam(sorted_bam_file_name, red_rem_bam_file_name, opt.redundancy, genome)
    print "Total reads: " + str(total) + "\nTotal retained reads: " + str(retained)

    # remove redundant reads in sorted control BAM file and count number of total reads and number of retained reads
    print "\nPreprocess the sorted control file to remove redundancy with threshold " + str(opt.redundancy) + "..."
    control_total, control_retained = SICER_MS.remove_redundant_reads_bam(sorted_control_file_name, red_rem_control_file_name, opt.redundancy, genome)
    print "Control file total reads: " + str(control_total) + "\nControl file total retained reads: " + str(control_retained) + "\n \n"

    # create HTSeq BAM readers that can iterate through all of the reads
    bam_iterator = HTSeq.BAM_Reader(red_rem_bam_file_name)
    control_bam_iterator = HTSeq.BAM_Reader(red_rem_control_file_name)

    print "Partition the genome in windows... \n"

    # evaluate if chip library is pair ended or single ended
    paired_end_bool = itertools.islice(bam_iterator,1).next().paired_end
    if paired_end_bool:
        # make dictionary of reads and windows and count total reads
        # read_dict: keys are chromosomes and values are a list of read positions
        # window_dict: keys are chromosomes and values are a list of window start coordinates for windows containing reads
        read_dict, window_dict, total_reads = SICER_MS.make_dict_of_reads_and_windows(bam_iterator, genome, opt.fragment_size,
                                                                             opt.window_size)
    elif not paired_end_bool:
        # make dictionary of reads and windows and count total reads
        # read_dict: keys are chromosomes and values are a list of read positions
        # window_dict: keys are chromosomes and values are a list of window start coordinates for windows containing reads
        read_dict, window_dict, total_reads = SICER_MS.make_dict_of_reads_and_windows(bam_iterator, genome, opt.fragment_size,
                                                                             opt.window_size)

    # make dictionary of reads and windows and count total reads for control file
    control_read_dict, control_window_dict, control_total_reads = SICER_MS.make_dict_of_reads_and_windows(control_bam_iterator,
                                                                                                 genome,
                                                                                                 opt.fragment_size,
                                                                                                 opt.window_size)

    print "Count reads in windows... \n"
    # get the read count and normalized read count of all windows in the bed file
    # create window counts dictionary window_counts_dict
    # add the window's score to the genomic array normalized_window_array
    # window_counts_dict: keys are chromosomes and values are a list of smaller
    # lists of the format [window_start, read_count, score] (the score will be calculated later)
    window_counts_dict, normalized_window_array = SICER_MS.get_window_counts_and_normalize(window_dict, read_dict, genome,
                                                                                  1000000, total_reads, opt.window_size)

    # get the read count and score of all windows in the control file file
    control_window_counts_dict, control_normalized_window_array = SICER_MS.get_window_counts_and_normalize(control_window_dict,
                                                                                                  control_read_dict,
                                                                                                  genome, 1000000,
                                                                                                  control_total_reads,
                                                                                                  opt.window_size)

    # write bedgraph file of normalized windows
    normalized_window_array.write_bedgraph_file(normalized_bedgraph_file_name)

    # write bedgraph file of normalized windows for control
    control_normalized_window_array.write_bedgraph_file(control_normalized_bedgraph_file_name)

    # write graph file for control reads
    SICER_MS.write_graph_file(control_window_counts_dict, opt.window_size, control_graph_file_name, genome)

    print "Find candidate islands exhibiting clustering... \n"

    # finds all islands using the dictionary of window counts and generates .scoreisland file
    # returns a genomic array island_array of all island tag counts and a list of islands (in dictionary format)
    # the dictionary keys of each island are 'island', 'score', and 'chip' (the read count)
    # also writes graph file
    island_array, islands_list = SICER_MS.find_islands(window_counts_dict, total_reads, opt.gap_size, opt.window_size, genome,
                                                opt.genome_fraction, e_value, score_island_file_name,
                                                graph_file_name, 2)


    # count the number of reads in the islands for both chip and control
    # returns updated list of islands including chip and control read counts and the total reads located in islands for
    # both island dictionaries
    islands_list, total_chip_reads_in_islands, total_control_reads_in_islands = SICER_MS.count_reads_in_islands(islands_list,
                                                                                                       read_dict,
                                                                                                       control_read_dict)

    print "Total chip reads in islands: " + str(total_chip_reads_in_islands)
    print "Total control reads in islands: " + str(total_control_reads_in_islands)

    # calculate the p-value and fold change (number of chip reads versus number of expected chip reads) for all islands
    # calculate alpha value for all islands
    # write island summary file
    # return list of islands islands_list; each island is a dictionary with keys 'island' (HTSeq genomic interval),
    # 'chip' (number of chip reads), 'control' (number of control reads), 'pvalue', 'fc' (fold change), and 'alpha'
    # also return HTSeq Genomic Array of all islands with their chip read count
    islands_list, island_array = SICER_MS.get_pvalue_fc_write_islandsummary(islands_list, total_reads, control_total_reads,
                                                     opt.genome_fraction, genome, islandsummary_file_name)

    print "\nIdentify significant islands using FDR criterion..."
    # given list of islands as dictionaries, filter all islands with alpha values meeting the significance threshold to write two files
    # write filtered island file (format: chr  start   end    chip_reads   control_reads   pvalue  fc  alpha)
    # write island bed file (format: chr   start   end     chip_reads)
    filtered_islands_list, filtered_island_array = SICER_MS.filter_islands_by_significance(islands_list,
                                                                                  filtered_island_file_name,
                                                                                  island_bed_file_name, FDR, genome)

    print "\nFilter reads with identified significant islands...\n"
    # given HTSeq bed_iterator and HTSeq Genomic Array that has chip read count assigned to all islands
    # finds all reads in the bed_iterator that are located in islands
    # if a read is located in an island, it is written to a bed file
    # creates a genomic array of all windows that have reads located in islands
    # returns a dictionary containing all reads located in islands and a dictionary containing all windows in islands
    # dictionary format: keys are chromosomes, values are sorted lists of all read/window positions
    islandfiltered_reads_dict, islandfiltered_windows_dict, total_chip_reads_in_islands = SICER_MS.filter_raw_tags_by_islands_bam(
                                                                                 bam_iterator,
                                                                                 filtered_island_array,
                                                                                 island_filtered_file_name,
                                                                                 opt.fragment_size, opt.window_size,
                                                                                 genome)

    # calculate the number of island filtered reads in all the windows comprising the islands
    # calculate normalized read count for each window
    # add the window's normalized read count to a genomic array (islandfilt_normalized_window_array)
    # the islandfilt_normalized_window_array will be used to write a bedgraph file
    islandfilt_window_counts_dict, islandfilt_normalized_window_array = SICER_MS.get_window_counts_and_normalize(islandfiltered_windows_dict,
                                                                                                        islandfiltered_reads_dict,
                                                                                                        genome,
                                                                                                        1000000,
                                                                                                        total_chip_reads_in_islands,
                                                                                                        opt.window_size)


    # write bedgraph file of normalized filtered islands
    islandfilt_normalized_window_array.write_bedgraph_file(islandfiltered_normalized_bedgraph_file_name)


if __name__ == "__main__":
    main(sys.argv)
