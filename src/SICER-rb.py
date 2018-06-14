# SICER-rb

# Version: 2.2
# Changes include modularization of essentials functions in SICER_MS script

import re, os, sys
from math import *
from string import *
from optparse import OptionParser
import time
import Background_island_probscore_statistics
import HTSeq
import bisect
import SICER_MS

def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--file", action="store", type="string", dest="file_name", metavar="<file>",
                      help="name of bed file (not including .bed extension)")
    parser.add_option("-g", "--genome", action="store", type="string", dest="genome_data", metavar="<file>",
                      help="name of reference genome (mm9 for mouse)")
    parser.add_option("-r", "--redundancy", action="store", type="int", dest="redundancy",
                      metavar="<file>", help="redundancy threshold")
    parser.add_option("-w", "--window_size", action="store", type="int", dest="window_size", metavar="<int>",
                      help="size of windows used to partition genome (200 for histones, 50 for TFs")
    parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>",
                      help="fragment size determines the shift (half of fragment_size of ChIP-seq read position, in bps)")
    parser.add_option("-p", "--genome_fraction", action="store", type="float", dest="genome_fraction", metavar="<int>",
                      help="effective genome fraction: 0.8 in most cases")
    parser.add_option("-s", "--gap_size", action="store", type="int", dest="gap_size", metavar="<int>",
                      help="maximum number of base pairs between windows in the same island (usually same as window size)")
    parser.add_option("-e", "--e-value", action="store", type="string", dest="e_value", metavar="<string>",
                      help="e-value used to determine significance")
    parser.add_option("-i", "--input_dir", action="store", type="string", dest="input_dir", metavar="<string>",
                      help="path to input directory")
    parser.add_option("-o", "--output_dir", action="store", type="string", dest="output_dir", metavar="<string>",
                      help="path to output directory")
    parser.add_option("-a", "--SICER_dir", action="store", type="string", dest="sicer_dir", metavar="<string>",
                      help="path to directory containing SICER files")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)

    #create string names for files
    #remove .bed extension
    file_name = opt.file_name[:-4]
    bed_file_name = opt.input_dir + "/" + opt.file_name
    sorted_bed_file_name = opt.output_dir + "/" + file_name + "_sorted_temp.bed"
    # This file stores the preprocessed raw bed file.
    red_rem_bed_file_name = opt.output_dir + "/" + file_name + "-" + str(opt.redundancy) + "-removed.bed"
    # This file stores the candidate islands.
    score_island_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + ".scoreisland"
    # This file stores the summary graph.
    graph_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + ".graph"
    # This file stores the island-filtered non-redundant raw reads
    island_filtered_file_name = opt.output_dir + "/" + file_name  + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + "-E" + str(opt.e_value) + "-islandfiltered.bed"
    # This file stores the sample summary graph in bedgraph format
    normalized_bedgraph_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-normalized.bedgraph"
    # This file stores normalized summary graph made by the island-filtered non-redundant raw reads in bedgraph format
    islandfiltered_normalized_bedgraph_file_name = opt.output_dir + "/" + file_name + "-W" + str(opt.window_size) + "-G" + str(opt.gap_size) + "-E" + str(opt.e_value) + "-islandfiltered-normalized.bedgraph"
    genome_file = opt.sicer_dir + "/genomes/" + opt.genome_data

    # read genome data from file containing genome data
    # store genome data in the dictionary genome
    genome = SICER_MS.get_genome_data(genome_file)

    # convert E_value to float
    e_value = float(opt.e_value)

    # sort bed file by chromosome, then by coordinate, then by strand
    print "\nSorting BED file..."
    SICER_MS.sort_bed_file(bed_file_name, sorted_bed_file_name)

    # remove redundant reads in bed file and count number of total reads and number of retained reads
    print "\nPreprocess the sorted BED file to remove redundancy with threshold " + str(opt.redundancy) + "..."
    total, retained = SICER_MS.remove_redundant_reads_bed(sorted_bed_file_name, red_rem_bed_file_name, opt.redundancy, genome)
    print "Total reads: " + str(total) + "\nTotal retained reads: " + str(retained) + "\n\n"

    # remove sorted bed file
    os.system('rm %s' % (sorted_bed_file_name))

    # create HTSeq bed reader that can iterate through all of the reads
    bed_iterator = HTSeq.BED_Reader(red_rem_bed_file_name)

    print "Partition the genome in windows... \n"
    # make dictionary of reads and windows and count total reads
    # read_dict: keys are chromosomes and values are a list of read positions
    # window_dict: keys are chromosomes and values are a list of window start coordinates for windows containing reads
    read_dict, window_dict, total_reads = SICER_MS.make_dict_of_reads_and_windows(bed_iterator, genome, opt.fragment_size,
                                                                         opt.window_size)

    print "Count reads in windows... \n"
    # calculate the number of island reads in all the windows comprising the islands
    # calculate normalized read count for each window
    # add the window's normalized read count to a genomic array (island_normalized_window_array)
    # the island_normalized_window_array will be used to write a bedgraph file
    window_counts_dict, normalized_window_array = SICER_MS.get_window_counts_and_normalize(window_dict, read_dict, genome,
                                                                                  1000000, total_reads, opt.window_size)
    # write bedgraph file of normalized islands
    normalized_window_array.write_bedgraph_file(normalized_bedgraph_file_name)

    print "Find candidate islands exhibiting clustering... \n"
    # finds all islands using the dictionary of window counts and generates .scoreisland file
    # returns a genomic array island_array of all island tag counts and a list of islands (in dictionary format)
    # the dictionary keys of each island are 'island', 'score', and 'chip' (the read count)
    # also writes graph file
    island_array, islands_list = SICER_MS.find_islands(window_counts_dict, total_reads, opt.gap_size, opt.window_size, genome,
                                                opt.genome_fraction, e_value, score_island_file_name,
                                                graph_file_name, 2)

    print "\nFilter reads with identified significant islands...\n"
    # given HTSeq bed_iterator and HTSeq Genomic Array that has chip read count assigned to all islands
    # finds all reads in the bed_iterator that are located in islands
    # if a read is located in an island, it is written to a bed file
    # creates a genomic array of all windows that have reads located in islands
    # returns a dictionary containing all reads located in islands and a dictionary containing all windows in islands
    # dictionary format: keys are chromosomes, values are sorted lists of all read/window positions
    islandfiltered_reads_dict, islandfiltered_windows_dict, total_reads_in_islands= SICER_MS.filter_raw_tags_by_islands(bed_iterator,
                                                                                 island_array,
                                                                                 island_filtered_file_name,
                                                                                 opt.fragment_size, opt.window_size,
                                                                                 genome)

    # calculate the number of island filtered reads in all the windows comprising the islands
    # calculate normalized read count for each window
    # add the window's normalized read count to a genomic array (islandfilt_normalized_window_array)
    # the islandfilt_normalized_window_array will be used to write a bedgraph file
    islandfiltered_window_counts_dict, islandfiltered_normalized_window_array = SICER_MS.get_window_counts_and_normalize(
                                                                                islandfiltered_windows_dict,
                                                                                islandfiltered_reads_dict, genome,
                                                                                1000000, total_reads_in_islands,
                                                                                opt.window_size)
    # write bedgraph file of normalized filtered islands
    islandfiltered_normalized_window_array.write_bedgraph_file(islandfiltered_normalized_bedgraph_file_name)


if __name__ == "__main__":
    main(sys.argv)
