# SICER_MS

# Version: 2.2
# Changes include modularization of essentials functions in SICER_MS script

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

########################## ALL SICER Code ##########################
########################## Code for get_genome_data ##########################
# module that reads genome data from a file with format: chromosome + /tab + length
# the module creates a dictionary of genome data and returns the dictionary
def get_genome_data(filename):
    genome_file = open(filename, 'r')
    genome_data = {}
    for line in genome_file:
        chrom = line.split()[0]
        length = int(line.split()[1])
        genome_data[chrom] = length
    genome_file.close()
    return genome_data
########################## End of Code for get_genome_data ##########################

########################## Code for make_dict_of_reads_and_windows ##########################
# make dictionary of reads and windows and count total reads
# read_dict: keys are chromosomes and values are a list of read positions
# window_dict: keys are chromosomes and values are a list of window start coordinates for windows containing reads
def make_dict_of_reads_and_windows(iterator, genome_data, fragment_size, window_size):

    total_reads = 0
    read_dict = {}
    window_dict = {}

    # add chromosomes to dictionaries as keys
    for chrom in genome_data:
        read_dict[chrom] = []
        window_dict[chrom] = []

    current_window_start = -1
    current_chrom = ""

    for read in iterator:
        #if not read.aligned:
        #    continue
        if read.iv.chrom in genome_data:
            # ensure that the read is not outside of the chromosome length
            if read.iv.start < 0:
                continue
            if read.iv.end >= genome_data[read.iv.chrom]:
                continue
            #calculate and store read position, if bed file or single end bam file use fragment_size
            read_pos = get_read_pos(read, fragment_size, genome_data)
						
            read_dict[read.iv.chrom].append(read_pos)
            total_reads += 1
            
            #round down to nearest window starting position
            window_start = read_pos / window_size * window_size
            window_end = window_start + window_size
            # make sure window is on chromosome
            if window_start >= 0 and window_end < genome_data[read.iv.chrom]:
                # make sure no duplicate windows are added
                #if window_start > current_window_start or current_chrom != read.iv.chrom:
                    window_dict[read.iv.chrom].append(window_start)
                    current_window_start = window_start
                    current_chrom = read.iv.chrom

    # sort read positions and window starts, and remove duplicate windows in the window dictionary
    for chrom in genome_data:

        # make sure read positions in read dictionary are sorted
        read_dict[chrom].sort()
        window_dict[chrom] = remove_duplicates_and_sort(window_dict[chrom])

    return read_dict, window_dict, total_reads
    

# remove duplicates in a list, sort the list, and return the sorted list with duplicates removed
def remove_duplicates_and_sort(l):
    no_dup = list(set(l))
    no_dup.sort()
    return no_dup

# determine read position given the fragment size
def get_read_pos(read, fragment_size, genome_data):
    shift = int(round(fragment_size / 2))
    if read.iv.strand == "+":
        read_pos = read.iv.start_d + shift
    elif read.iv.strand == "-":
        read_pos = read.iv.start_d - shift
	# make sure read_pos is located on chromosome
    if read_pos < 0:
        read_pos = 0
    if read_pos >= genome_data[read.iv.chrom]:
        read_pos = genome_data[read.iv.chrom] - 1
    return read_pos



########################## End of Code for make_dict_of_reads_and_windows ##########################

########################## Code for get_window_counts_and_normalize ##########################
# calculate the read count of all windows in a dictionary of windows using the dictionary of read tags
# create window counts dictionary (keys are chromosomes and values are a list of lists where each small
# list is of the format [window start, read count, score]; the score will be calculated in another module
# add the window's normalized read count to the genomic array normalized_window_array
# return the window counts dictionary and the genomic array normalized_window_array
def get_window_counts_and_normalize(window_dict, tags_dict, genome_data, scaling_factor, total_reads, window_size):
    # dictionary to store read count in each window
    window_counts_dict = {}
    # HTSeq genomic array to store normalized score for each window (used to generate bedgraph file)
    normalized_window_array = HTSeq.GenomicArray(genome_data, stranded=False, typecode='d')

    # create chromosome keys in window counts dictionary for all chromosomes in genome; the values are empty lists
    for chrom in genome_data:
        window_counts_dict[chrom] = []

    # iterate through all chromosomes in the genome
    for chrom in genome_data:
        # iterate through all windows on the chromosome
        for window_start in window_dict[chrom]:
            # get read count in window
            read_count = get_read_count_in_window(chrom, window_start, window_size, tags_dict)

            window_counts_dict[chrom].append([window_start, read_count, 0])

            # calculate normalized read count
            normalized_count = float(read_count) * float(scaling_factor) / float(total_reads)
            window_end = window_start + window_size
            window = HTSeq.GenomicInterval(chrom, window_start, window_end)
            # assign normalized read count to window on HTSeq genomic array
            normalized_window_array[window] = normalized_count

    return window_counts_dict, normalized_window_array

# module that calculates the read count in a window
# requires the chromosome, window start, window size, and a dictionary of sorted lists of read positions
# returns the read count
def get_read_count_in_window(chrom, window_start, window_size, tags_dict):
    # subtract one from the window end because the windows are half open
    window_end = window_start + window_size - 1
    # get the index of the first read in the window
    start_ind = bisect.bisect_left(tags_dict[chrom], window_start)
    # get the index of the first read outside of the window
    end_ind = bisect.bisect_right(tags_dict[chrom], window_end)
    # the read count is equal to the index of the first read outside the window minus the index of the first read in the window
    read_count = end_ind - start_ind
    return read_count
########################## End of Code for get_window_counts_and_normalize ##########################

########################## Code for find_islands ##########################
# finds all islands using the dictionary of window counts and generates .scoreisland file
# returns a genomic array island_array of all island tag counts and a list of islands (in dictionary format)
# the dictionary keys of each island are 'island', 'score', and 'chip' (the read count)
# also writes a graph file
def find_islands(window_counts_dict, total_reads, gap, window_size, genome_data,
                 fraction, evalue, scoreisland_file_name, graph_file_name, window_size_buffer = 2):

    # calculate length of genome (sum up lengths of all chromosomes)
    genome_length = sum(genome_data.values())
    print "Genome length: " + str(genome_length)

    # calculate effective genome length by multiplying by the effective genome fraction
    genome_length = int(fraction * genome_length)
    print "Effective Genome length: " + str(genome_length)

    # calculate average number of reads in a window
    average = float(total_reads) * window_size / genome_length
    print "Window average: " + str(average)

    window_pvalue = 0.20
    print "Window p-value: " + str(window_pvalue)

    bin_size = 0.001
    background = Background_island_probscore_statistics.Background_island_probscore_statistics(total_reads,
                                                                                               window_size, gap,
                                                                                               window_pvalue,
                                                                                               genome_length,
                                                                                               bin_size)
    # determine minimum number of tags required in a qualified window
    min_tags_in_window = background.min_tags_in_window
    print "Minimum number of tags in a qualified window: " + str(min_tags_in_window)

    # determine score threshold for an island
    score_threshold = background.find_island_threshold(evalue)
    print "Score threshold from random background: " + str(score_threshold)

    print "Mark and write islands..."

    # iterate through all windows and create graph file
    # see if each window has minimum number of tags required and calculate each window's score
    # windows with minimum number of tags are retained in filtered_windows_dict
    # filtered_windows_dict: keys are chromosomes and values are a list of smaller lists of the format [window_start, read_count, score]
    filtered_windows_dict = write_graph_file_and_filter_windows(window_counts_dict, min_tags_in_window, window_size,
                                                           graph_file_name, average, genome_data)

    # generate list of islands that contains the island, its read count, and its score
    islands_list = form_islands(filtered_windows_dict, window_size, gap, window_size_buffer,
                                score_threshold, genome_data)

    # genomic array to store final islands (value = read count if island present)
    island_array = HTSeq.GenomicArray(genome_data, stranded=False, typecode='d')

    score_island_file = open(scoreisland_file_name, 'w')

    # write all final islands to a bed file (.scoreisland file) and add them to island_array
    for i in islands_list:
        # line for .score_island file
        outline = i['island'].chrom + "\t" + str(i['island'].start) + "\t" + str(i['island'].end) + \
                  "\t" + str(i['score']) + "\n"
        score_island_file.write(outline)

        # assign island's read count to the island on the HTSeq Genomic Array
        island_array[i['island']] = i['chip']

    score_island_file.close()

    print "Total number of islands: " + str(len(islands_list)) + "\n"

    return island_array, islands_list

# create a graph file of all windows and their read count using a window counts dictionary
# window_counts_dict: keys are chromosomes and values are a list of smaller lists of the format [window_start, read_count, score] (score may not be calculated yet)
# retain only those windows that satisfy minimum number of tags in the dictionary filtered_windows_dict and calculate their score
# filtered_windows_dict: keys are chromosomes and values are a list of smaller lists of the format [window_start, read_count, score]
# return filtered_windows_dict
def write_graph_file_and_filter_windows(window_counts_dict, min_tags_in_window, window_size, graph_file_name,
                                        average, genome_data):

    # dictionary that will contain all windows satisfying minimum number of tags
    filtered_windows_dict = {}
    for chrom in genome_data:
        filtered_windows_dict[chrom] = []

    graph_file = open(graph_file_name, 'w')

    # iterate through all chromosomes in the genome
    for chrom in genome_data:
        # make sure chromosome has windows in it
        if len(window_counts_dict[chrom]) > 0:
            # iterate through all windows on the chromosome
            for window in window_counts_dict[chrom]:
                read_count = window[1]
                window_start = window[0]
                window_end = window_start + window_size
                # line to add to graph file
                graph_file_line = chrom + '\t' + str(window_start) + '\t' \
                                  + str(window_end) + '\t' + str(int(read_count)) + '\n'
                graph_file.write(graph_file_line)

                # see if window has minimum number of reads
                if read_count < min_tags_in_window:
                    # if read count doesn't meet threshold, window score is -1
                    score = -1
                else:
                    # calculate probability of window using poisson function
                    prob = poisson(read_count, average)

                    # calculate score of window based on the probability
                    if prob < 1e-250:
                        score = 1000
                    else:
                        score = -log(prob)

                # assign the score to the window's information list of format [window start, read count, score]
                window[2] = score

                # if window contains minimum number of tags required, append to filtered_windows list along with
                # its calculated score
                if score > 0:
                    filtered_windows_dict[chrom].append(window)

    graph_file.close()
    return filtered_windows_dict

def fact(m):
    value = 1.0
    if m != 0:
        while m != 1:
            value = value * m
            m = m - 1
    return value

def factln(m):
    if m < 20:
        value = 1.0
        if m != 0:
            while m != 1:
                value = value*m
                m = m - 1
        return log(value)
    else:
        return m*log(m) - m + log(m*(1+4*m*(1+2*m)))/6.0 + log(pi)/2


def poisson(i, average):
    if i < 20:
        return exp(-average) * average**i / fact(i)
    else:
        exponent = -average + i * log(average) - factln(i)
        return exp(exponent)

# given a list of windows, gap size, window size buffer, and score threshold,
# return a list of the islands that meet the score threshold
# the returned list contains a list of all the islands in dictionary format
# the dictionary keys of each island are 'island', 'score', and 'chip' (the read count)
def form_islands(windows_dict, window_size, gap, window_size_buffer, score_threshold, genome_data):

    # maximum distance between windows in same island
    proximal_island_dist = gap + window_size_buffer

    # list to store all potential islands (that may or may not meet score threshold)
    islands = []

    # iterate through all the chromosomes in the genome
    for chrom in genome_data:

        # make sure the chromosome has windows in it
        if len(windows_dict[chrom]) > 0:

            # get the information of the first window on the chromosome
            current_island_start = windows_dict[chrom][0][0]
            current_island_end = current_island_start + window_size
            # create HTSeq Genomic Interval object to store island information
            current_island_interval = HTSeq.GenomicInterval(chrom, current_island_start, current_island_end)
            # island has format [HTSeqGenomicInterval, read count, score]
            current_island = [current_island_interval, windows_dict[chrom][0][1], windows_dict[chrom][0][2]]

            # iterate through all windows and join windows together to form islands
            for window in windows_dict[chrom][1:]:
                window_start = window[0]
                window_end = window_start + window_size

                # calculate distance between window start and end of previous island
                dist = window_start - current_island[0].end

                # if distance is less than proximal island distance
                if dist <= proximal_island_dist:
                    # change end of island to be the end of the window
                    new_interval = HTSeq.GenomicInterval(chrom, current_island[0].start, window_end)
                    current_island[0] = new_interval

                    # add window's read count and score to the island
                    current_island[1] += window[1]
                    current_island[2] += window[2]

                # if distance is greater than proximal island distance
                else:
                    # append the current island to the list of islands
                    islands.append(current_island)

                    # the window will be the new island
                    new_island_interval = HTSeq.GenomicInterval(chrom, window_start, window_end)
                    current_island = [new_island_interval, window[1], window[2]]

            # append the last island
            islands.append(current_island)

    # list to store all islands that meet score threshold
    final_islands = []

    # iterate through all islands and retain only those islands that meet score threshold
    for island in islands:
        if island[2] >= score_threshold - 0.0000000001:
            # create a dictionary object for the island
            # 'island' is the HTSeq.GenomicInterval, 'score' is the island's score, and 'chip' is the read count
            island_dict = {'island': island[0], 'score': island[2], 'chip': island[1]}

            # append island dictionary to the list of final islands
            final_islands.append(island_dict)

    return final_islands
########################## End of Code for find_islands ##########################
########################## End of ALL SICER Code ##########################

########################## BED SICER Code ##########################
########################## Code for sort_bed_file ##########################
# module that sorts a BED file first by chromosome, then by coordinate, and then by strand
def sort_bed_file(infile, outfile):

    f = open(infile, 'r')
    o = open(outfile, 'w')

    try:
        # sort first by chromosome, then by coordinate, then by strand
        if os.system('sort -k1,1 -k2,3n -k6,6 %s > %s' % (infile, outfile)):
            raise
    except:
        sys.stderr.write("reads do not exist in " + str(infile) + "\n")
########################## End of Code for sort_bed_file ##########################

########################## Code for remove_redundant_reads_bed ##########################
# module that removes redundant lines in a BED file
# BED file must be sorted first
def remove_redundant_reads_bed(infile, outfile, cutoff, genome_data):

    f = open(infile, 'r')
    o = open(outfile, 'w')

    current_start = 0
    current_end = 0
    current_count = 1
    current_strand = "!"
    # count total number of reads
    total = 0
    # count number of retained reads (number of reads with duplicates removed)
    retained = 0

    for line in f:
        line = line.strip()
        # make sure line is an actual read and not comments
        if not re.match("#", line):

            total += 1

            chr = str(line.split()[0])
            # make sure chr is in genome
            if chr not in genome_data:
                continue

            start = int(line.split()[1])

            # make sure read is located on genome
            if start < 0:
                print "Illegitimate read with start less than zero is ignored: " + line
                continue
            end = int(line.split()[2])
            if end > genome_data[chr]:
                print "Illegitimate read with end beyond chromosome length " + str(genome_data[chr]) + " is ignored: " \
                      + line
                continue

            strand = str(line.split()[5])

            # see if read is a duplicate
            if start != current_start:
                o.write('\t'.join(line.split())+'\n')
                retained += 1
                current_start = start
                current_end = end
                current_strand = strand
                current_count = 1
            elif end != current_end:
                o.write('\t'.join(line.split())+'\n')
                retained += 1
                current_start = start
                current_end = end
                current_strand = strand
                current_count = 1
            elif strand != current_strand:
                o.write('\t'.join(line.split())+'\n')
                retained += 1
                current_start = start
                current_end = end
                current_strand = strand
                current_count = 1

            # read is a duplicate
            else:
                current_count += 1
                assert current_start == start
                assert current_end == end
                assert current_strand == strand

                if current_count <= cutoff:
                    o.write('\t'.join(line.split())+'\n')
                    retained += 1

    f.close()
    o.close()

    # return total number of reads and total number of retained reads
    return total, retained
########################## End of Code for remove_redundant_reads_bed ##########################

########################## Code for filter_raw_tags_by_islands ##########################
# given HTSeq bed_iterator and HTSeq Genomic Array that has chip read count assigned to all islands
# finds all reads in the bed_iterator that are located in islands
# if a read is located in an island, it is written to a bed file
# creates a genomic array of all windows that have reads located in islands
# returns a dictionary containing all reads located in islands and a dictionary containing all windows in islands
# dictionary format: keys are chromosomes, values are sorted lists of all read/window positions
def filter_raw_tags_by_islands(bed_iterator, island_array, filename, fragment_size, window_size, genome_data):

    file = open(filename, 'w')

    # count total number of reads in islands
    total_reads_in_islands = 0
    window_counts = HTSeq.GenomicArray(genome_data, stranded=False, typecode='d')

    # dictionary to store reads located in islands
    islandfiltered_reads_dict = {}
    # dictionary to store windows located in islands
    islandfiltered_windows_dict = {}

    for chrom in genome_data:
        islandfiltered_reads_dict[chrom] = []
        islandfiltered_windows_dict[chrom] = []

    # iterate through all reads in the bed iterator and write to file if read lands in an island
    for read in bed_iterator:

        # apply shift to determine coordinate
        read_pos = get_read_pos(read, fragment_size, genome_data)
        position = HTSeq.GenomicPosition(read.iv.chrom, read_pos)

        # determine if read lands in an island by using the genomic array island_array; if it does, write to file
        if island_array[position] >= 1:
            # line to add to BED file
            line = str(read.iv.chrom) + "\t" + str(read.iv.start) + "\t" + str(read.iv.end) \
                    + "\t" + str(read.name) + "\t" + str(read.score) + "\t" + str(read.iv.strand) + "\n"
            file.write(line)
            total_reads_in_islands += 1
            islandfiltered_reads_dict[read.iv.chrom].append(read_pos)
            # determine window position of read
            window_start = read_pos / window_size * window_size
            if window_start < 0:
                window_start = 0
            islandfiltered_windows_dict[read.iv.chrom].append(window_start)

    file.close()

    # sort and remove duplicates for reads dictionary and windows dictionary
    for chrom in genome_data:
        islandfiltered_reads_dict[chrom].sort()
        islandfiltered_windows_dict[chrom] = remove_duplicates_and_sort(islandfiltered_windows_dict[chrom])

    return islandfiltered_reads_dict, islandfiltered_windows_dict, total_reads_in_islands
########################## End of Code for filter_raw_tags_by_islands ##########################
########################## End of BED SICER Code ##########################

########################## BAM Paired to Single Read Code ##########################
def write_to_outfile(outfile, GI, read):
    outfile.write(str(GI.chrom)+"\t"+str(GI.start)+"\t"+str(GI.end)+"\t"+str(read.name.split('/')[0])+"\t"+str(read.score)+"\t"+str(GI.strand)+"\n")

def paired_to_single_read(bamfile):

    #Sort paired-end BAM file by name so that paired reads are adjacent
    bam_reader = HTSeq.BAM_Reader(bamfile)
    os.system('samtools sort -O BAM -n %s > %s' % (bamfile, bamfile[:-4]+"_sorted.bam"))

    #BAM file is converted to BED format
    os.system('bamToBed -i %s > %s' % (bamfile[:-4]+"_sorted.bam", bamfile[:-4]+"_sorted.bed"))

    #BED iterator created to traverse BED file
    bed_iterator = HTSeq.BED_Reader(bamfile[:-4]+"_sorted.bed")

    outfile = open(bamfile[:-4]+"_single.bed", 'w')

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
                write_to_outfile(outfile, singleRead_iv, read)
                finalcount = finalcount + 1
            else:
                #print "Error: paired reads not on same chromosome."
                error_count = error_count + 1



    #print "pairedRead to singleRead conversion is done running.\nThere were " + str(error_count) + " pairs on different chromosomes"
    #print "Reads written to outfile: " + str(finalcount)
########################## End BAM Paired to Single Read Code ##########################



########################## BAM SICER Code ##########################
########################## Code for remove_redundant_reads_bam ##########################
# Removes duplicates and unmapped reads in a singled ended bam file
# Bam file passed to this module is presorted
# The module traverses bam file containing reads using bam_reader
# Module first checks to see if the read is unmapped or belongs to chrM and filters it out if it is
# Module then checks if read is already present in output .bam file within the boundaries of the specified redundancy threshold and filters it out if necessary
def remove_redundant_reads_bam(file_name, out_file_name, cutoff, genome_data):
    bam_reader = HTSeq.BAM_Reader(file_name)
    bam_writer = HTSeq.BAM_Writer.from_BAM_Reader(out_file_name, bam_reader)

    total = 0
    retained = 0
    current_read = None
    earlyDup = None

    trashDup = 0
    onlyTrash = 0
    removedDup = 0


    for read in bam_reader:
        total += 1
        # checking first read
        if current_read is None:
            if isTrash(read, genome_data) == False:
                if total == 1:
                    print "First read is normal"
                current_read = read
                bam_writer.write(read)
                retained += 1
                current_count = 1
                continue
            else:
                if earlyDup is None:
                    onlyTrash+=1
                    earlyDup = read
                    continue
                elif areDuplicates(read, earlyDup):
                    if trashDup == 0:
                        onlyTrash = 0
                        trashDup = 2
                    else:
                        trashDup+=1
                    earlyDup = read
                    continue
                else:
                    onlyTrash+=1
                    earlyDup = read
                    continue

        # check if read is aligned or reads of chrM or read not in genome or reads are duplciates
        # aka read is trash
        if isTrash(read, genome_data):
            onlyTrash+=1
            continue

        # read is not a duplicate
        if areDuplicates(read, current_read) == False:
            bam_writer.write(read)
            retained += 1
            current_read = read
            current_count = 1
            continue

        # read is a duplicate of current_read
        else:
            current_count += 1
            current_read = read

            if current_count <= cutoff:
                bam_writer.write(read)
                retained += 1
            else:
                removedDup+=1

    bam_writer.close()

    print "Complete number of trash reads is: " + str(onlyTrash)
    print "Total number of removed non-trash redundancies is: " + str(removedDup)

    return total, retained

def isTrash(read, genome_data):
    if ((read.aligned == False) or (read.iv.chrom not in genome_data) or (read.iv.start < 0) or (read.iv.end >= genome_data[read.iv.chrom])):
        return True
    else:
        return False

def areDuplicates(readone, readtwo):
    if ((readone.iv.start == readtwo.iv.start) and (readone.iv.end == readtwo.iv.end) and (readone.iv.strand == readtwo.iv.strand)):
        return True
    else:
        return False
########################## End of Code for remove_redundant_reads_bam ##########################

########################## Code for filter_raw_tags_by_islands_bam ##########################
# given HTSeq bed_iterator and HTSeq Genomic Array that has chip read count assigned to all islands
# finds all reads in the bed_iterator that are located in islands
# if a read is located in an island, it is written to a bed file
# creates a genomic array of all windows that have reads located in islands
# returns a dictionary containing all reads located in islands and a dictionary containing all windows in islands
# dictionary format: keys are chromosomes, values are sorted lists of all read/window positions
def filter_raw_tags_by_islands_bam(iterator, island_array, filename, fragment_size, window_size, genome_data):

    file = open(filename, 'w')

    # count total number of reads in islands
    total_reads_in_islands = 0
    window_counts = HTSeq.GenomicArray(genome_data, stranded=False, typecode='d')

    # dictionary to store reads located in islands
    islandfiltered_reads_dict = {}
    # dictionary to store windows located in islands
    islandfiltered_windows_dict = {}

    for chrom in genome_data:
        islandfiltered_reads_dict[chrom] = []
        islandfiltered_windows_dict[chrom] = []

    # iterate through all reads in the bam iterator and write to file if read lands in an island
    for read in iterator:

        if not read.aligned:
            continue
        if read.iv.chrom not in genome_data:
            continue
        # apply shift to determine coordinate
        read_pos = get_read_pos(read, fragment_size, genome_data)
        position = HTSeq.GenomicPosition(read.iv.chrom, read_pos)

        if position.pos >= genome_data[read.iv.chrom]:
            position.pos = genome_data[read.iv.chrom] - 1

        # determine if read lands in an island by using the genomic array island_array; if it does, write to file
        if island_array[position] >= 1:
            # line to add to BED file
            line = str(read.iv.chrom) + "\t" + str(read.iv.start) + "\t" + str(read.iv.end) \
                    + "\t" + str(read.read.name) + "\t" + str(read.aQual) + "\t" + str(read.iv.strand) + "\n"
            file.write(line)
            total_reads_in_islands += 1
            islandfiltered_reads_dict[read.iv.chrom].append(read_pos)
            # determine window position of read
            window_start = read_pos / window_size * window_size
            if window_start < 0:
                window_start = 0
            islandfiltered_windows_dict[read.iv.chrom].append(window_start)

    file.close()

    # sort and remove duplicates for reads dictionary and windows dictionary
    for chrom in genome_data:
        islandfiltered_reads_dict[chrom].sort()
        islandfiltered_windows_dict[chrom] = remove_duplicates_and_sort(islandfiltered_windows_dict[chrom])

    return islandfiltered_reads_dict, islandfiltered_windows_dict, total_reads_in_islands
########################## End of Code for filter_raw_tags_by_islands_bam ##########################
########################## End of BAM SICER Code ##########################

########################## NON-RB SICER Code ##########################
########################## Code for write_graph_file ##########################
# create a graph file given a window counts dictionary, the window size, the name of the file, and the genome data
# window_counts_dict: keys are chromosomes and values are a list of smaller lists of the format [window_start, read_count]
def write_graph_file(window_counts_dict, window_size, graph_file_name, genome_data):

    graph_file = open(graph_file_name, 'w')

    # iterate through all chromosomes in genome
    for chrom in genome_data:

        # make sure chromosome has windows in it
        if len(window_counts_dict[chrom]) > 0:

            # iterate through all windows on chromosome
            for window in window_counts_dict[chrom]:

                read_count = window[1]
                window_start = window[0]
                window_end = window_start + window_size

                # line to write to graph file
                graph_file_line = chrom + '\t' + str(window_start) + '\t' \
                                  + str(window_end) + '\t' + str(int(read_count)) + '\n'

                graph_file.write(graph_file_line)

    graph_file.close()
########################## End of Code for write_graph_file ##########################

########################## Code for count_reads_in_islands ##########################
# given a list of islands, and two dictionaries of read positions (for chip and control, where keys are chromosomes and values are sorted lists of read positions)
# counts the number of reads in the islands from both dictionaries
# returns updated list of islands including read counts and the total reads located in islands for both island dictionaries
def count_reads_in_islands(islands_list, chip_tags_dict, control_tags_dict):

    total_chip_reads_in_islands = 0
    total_control_reads_in_islands = 0

    # iterate through all islands in list of islands
    for island in islands_list:

        # calculate chip read count in island
        chip_count = get_read_count_in_region(island['island'], chip_tags_dict)
        # calculate control read count in island
        control_count = get_read_count_in_region(island['island'], control_tags_dict)

        # create dictionary keys for chip and control read counts
        island['chip'] = chip_count
        island['control'] = control_count

        # add chip and control read counts to total
        total_chip_reads_in_islands += chip_count
        total_control_reads_in_islands += control_count

    return islands_list, total_chip_reads_in_islands, total_control_reads_in_islands
########################## End of Code for count_reads_in_islands ##########################

########################## Code for get_read_count_in_region ##########################
# module that calculates the readcount in an HTSeq genomic interval object (region)
# requires the region and a dictionary of sorted lists of read positions (the keys are the chromosomes)
# returns the read count in the region
def get_read_count_in_region(region, tags_dict):

    # need to subtract 1 from region end coordinate because we are using half open notation
    region_end = region.end - 1
    # get the index of the first read in the region
    start_ind = bisect.bisect_left(tags_dict[region.chrom], region.start)
    # get the index of the first read outside of the region
    end_ind = bisect.bisect_right(tags_dict[region.chrom], region_end)
    # the read count is equal to the index of the first read outside of the region minus the index of the first read inside the region
    read_count = end_ind - start_ind

    return read_count
########################## End of Code for get_read_count_in_region ##########################


########################## Code for get_pvalue_fc_write_islandsummary ##########################
# calculate the p-value and fold change (number of chip reads versus number of expected chip reads) for all islands
# calculate alpha value for all islands
# write island summary file
# return list of islands; each island is a dicitonary with keys 'island' (HTSeq genomic interval), 'chip' (number of
# chip reads), 'control' (number of control reads), 'pvalue', 'fc' (fold change), and 'alpha'
# also return HTSeq Genomic Array of all islands with their chip read count
def get_pvalue_fc_write_islandsummary(islands_list, chip_library_size, control_library_size, genome_fraction,
                                      genome_data, filename):

    # list to store p-values
    pvalue_list = []

    # list to store island information
    result_list = []

    #HTSeq Genomic Array to store all of the chip read counts in each island
    island_array = HTSeq.GenomicArray(genome_data, stranded=False, typecode='d')

    # calculate scaling factor
    scaling_factor = chip_library_size * 1.0 / control_library_size

    # calculate genome size by summing all chromosome lengths
    genome_size = sum(genome_data.values())

    # change genome size to be the effective genome size by multiplying by effective genome fraction
    genome_size *= genome_fraction

    # iterate through all of the islands
    for island in islands_list:
        # get number of chip reads in island
        observation = island['chip']
        # get number of control reads in island
        control_tag = island['control']

        if island['control'] > 0:
            # calculate expected number of reads (average) in the island
            average = control_tag * scaling_factor
            # calculate the fold change (number of chip reads divided by average expected reads)
            fc = float(observation) / float(average)

        # if there are no control reads in the island
        else:
            # calculate the length in base pairs of the island
            length = island['island'].end - island['island'].start
            # estimate the average expected reads in the island by calculating the average number of control reads in a region of the length of the island
            average = length * control_library_size * 1.0 / genome_size
            # multiply the average by the scaling factor (but if the average is less than 0.25, just use 0.25)
            average = min(0.25, average) * scaling_factor
            # calculate the fold change (number of chip reads divided by average expected reads)
            fc = float(observation) / float(average)

        # calculate p-value given the number chip reads and the average expected reads
        if observation > average:
            pvalue = scipy.stats.poisson.sf(island['chip'], average)[()]
        # if number of chip reads is not greater than average expected reads, set pvalue to 1
        else:
            pvalue = 1

        # append p-value to pvalue_list
        pvalue_list.append(pvalue)

        # dictionary that will contain all island information
        item_dic = {}
        item_dic['chrom'] = island['island'].chrom
        item_dic['start'] = island['island'].start
        item_dic['end'] = island['island'].end
        item_dic['chip'] = observation
        item_dic['control'] = control_tag
        item_dic['pvalue'] = pvalue
        item_dic['fc'] = fc

        # append dicitonary of island information (item_dic) to list result_list
        result_list.append(item_dic)

        # assign the chip read count to the island region on the genomic array
        island_array[island['island']] = island['chip']

    # list to store islands once alpha is calculated
    new_islands_list = []
    islandsummary_file = open(filename, 'w')
    # create scipy array of p-values
    pvaluearray = scipy.array(pvalue_list)
    # rank the p-values using scipy
    pvaluerankarray = scipy.stats.rankdata(pvaluearray)
    totalnumber = len(result_list)

    # iterate through the islands in the result_list list to calculate alpha and write island summary file
    for i in range(totalnumber):
        item = result_list[i]
        # calculate alpha value
        alpha = pvalue_list[i] * totalnumber / pvaluerankarray[i]
        if alpha > 1:
            alpha = 1

        # line to write to island summary file
        outline = item['chrom'] + "\t" + str(item['start']) + "\t" + str(item['end']) + "\t" + str(
            item['chip']) + "\t" + str(item['control']) + "\t" + str(item['pvalue']) + "\t" + str(
            item['fc']) + "\t" + str(alpha) + "\n"
        # write line to island summary file
        islandsummary_file.write(outline)
        interval = HTSeq.GenomicInterval(item['chrom'], item['start'], item['end'])
        # new dictionary for island (now including the calculated alpha value)
        updated_island = {'island': interval, 'chip': item['chip'], 'control': item['control'],
                          'pvalue': item['pvalue'], 'fc': item['fc'], 'alpha': alpha}

        # append island dictionary to new_islands_list (now includes alpha calculation)
        new_islands_list.append(updated_island)

    islandsummary_file.close()

    return new_islands_list, island_array
########################## End of Code for get_pvalue_fc_write_islandsummary ##########################

########################## Code for filter_islands_by_significance ##########################
# given list of islands as dictionaries, filter all islands with alpha values meeting the significance threshold to write two files
# write filtered island file (format: chr  start   end    chip_reads   control_reads   pvalue  fc  alpha)
# write island bed file (format: chr   start   end     chip_reads)
def filter_islands_by_significance(islands_list, filtered_island_filename, island_bed_file_name, significance,
                                   genome_data):

    filtered_island_file = open(filtered_island_filename, 'w')
    filtered_bed_file = open(island_bed_file_name, 'w')

    filtered_island_array = HTSeq.GenomicArray(genome_data, stranded=False, typecode='d')

    # list to store significant islands
    filtered_island_list = []

    # count number of significant islands
    totalislands = 0

    # iterate through all islands
    for island in islands_list:

        # only use islands who's alpha value is less than significance
        if island['alpha'] <= significance:

            totalislands += 1

            filtered_island_line = island['island'].chrom + "\t" + str(island['island'].start) + "\t" + \
                                   str(island['island'].end) + "\t" + str(island['chip']) + "\t" + \
                                   str(island['control']) + "\t" + str(island['pvalue']) + "\t" + str(island['fc']) + \
                                   "\t" + str(island['alpha']) + "\n"

            filtered_island_bed_line = island['island'].chrom + "\t" + str(island['island'].start) + "\t" + \
                                       str(island['island'].end) + "\t" + str(island['chip']) + "\n"

            filtered_island_file.write(filtered_island_line)
            filtered_bed_file.write(filtered_island_bed_line)

            filtered_island_list.append(island)
            filtered_island_array[island['island']] = island['chip']


    print "Given significance " + str(significance) + ", there are " + str(totalislands) + " significant islands"

    filtered_island_file.close()
    filtered_bed_file.close()

    return filtered_island_list, filtered_island_array
########################## End of Code for filter_islands_by_significance ##########################
########################## End of NON-RB SICER Code ##########################
