import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator

import HTSeq
from SICER import make_dict_of_reads
from SICER import get_genome_data
from SICER import get_read_count_in_region

import scipy.stats


def pvalue(chip_read_count, control_read_count, scaling_factor, pseudo_count):
    """
    Currently using poisson distribution

    scaling_factor: the factor that accounts for the differences of control library and ChIP library. effective control read count
    is control_read_count * scaling factor
    pseudocount: when control_read_count is zero, replace zero with pseudocount to alleviate the impact of statistical fluctuation

    output: pvalue
    """
    if control_read_count > 0:
        average = control_read_count * scaling_factor
    else:
        average = pseudo_count * scaling_factor
    if chip_read_count > average:
        pvalue = scipy.stats.poisson.sf(chip_read_count, average)[()]
    else:
        pvalue = 1
    return pvalue


def fdr(pvalue_list):
    """
    Calculate the multiple testing corrected p-value using BH
    """
    fdr_list = []
    pvaluearray = scipy.array(pvalue_list)
    totalnumber = len(pvalue_list)
    pvaluerankarray = scipy.stats.rankdata(pvaluearray)
    for i in range(totalnumber):
        fdr_value = pvalue_list[i] * totalnumber / pvaluerankarray[i]
        if fdr_value > 1:
            fdr_value = 1
        fdr_list.append(fdr_value)
    return fdr_list


def main(argv):
    parser = OptionParser()
    parser.add_option("-g", "--genome", action="store", type="string", dest="genome_data", help="species, mm9, hg18, etc",
                      metavar="<str>")
    parser.add_option("-a", "--rawreadfileA", action="store", type="string", dest="readfileA", metavar="<file>",
                      help="raw read file A in bed format")
    parser.add_option("-b", "--rawreadfileB", action="store", type="string", dest="readfileB", metavar="<file>",
                      help="raw read file B in bed format")
    parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>",
                      help="average size of a fragment after A experiment")
    parser.add_option("-d", "--islandfile", action="store", type="string", dest="islandfile", metavar="<file>",
                      help="island file in BED format")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>",
                      help="island read count summary file")


    (opt, args) = parser.parse_args(argv)
    if len(argv) < 12:
        parser.print_help()
        sys.exit(1)

    # create HTSeq BED_Readers for BED files
    file_A_iterator = HTSeq.BED_Reader(opt.readfileA)
    file_B_iterator = HTSeq.BED_Reader(opt.readfileB)
    island_file_iterator = HTSeq.BED_Reader(opt.islandfile)

    genome_file = opt.sicer_dir + "/genomes/" + opt.genome_data
    genome = get_genome_data(genome_file)

    read_dict_A, A_library_size = make_dict_of_reads(file_A_iterator, genome, opt.fragment_size)

    read_dict_B, B_library_size = make_dict_of_reads(file_B_iterator, genome, opt.fragment_size)

    print "Library size of " + opt.readfileA + ":  " + str(A_library_size)
    print "Library size of " + opt.readfileB + ":  " + str(B_library_size)

    A_reads_in_islands = 0
    B_reads_in_islands = 0

    islands_list = []

    island_A_readcount_list = []
    island_B_readcount_list = []

    # Find read counts on the islands
    for region in island_file_iterator:

        read_count_A = get_read_count_in_region(region.iv, read_dict_A)
        A_reads_in_islands += read_count_A
        island_A_readcount_list.append(read_count_A)
        read_count_B = get_read_count_in_region(region.iv, read_dict_B)
        B_reads_in_islands += read_count_B
        island_B_readcount_list.append(read_count_B)

        island = {'region': region.iv, 'A_count': read_count_A, 'B_count': read_count_B}
        islands_list.append(island)

    pvalue_A_vs_B_list = []
    pvalue_B_vs_A_list = []


    print "Total number of A reads on islands is: " + str(A_reads_in_islands)
    print "Total number of B reads on islands is: " + str(B_reads_in_islands)

    library_scaling_factor = A_library_size * 1.0 / B_library_size
    pseudo_count = 1
    pvalue_A_vs_B_list = []
    pvalue_B_vs_A_list = []

    # Calculate the p value.
    for island in islands_list:
        A_count = island['A_count']
        B_count = island['B_count']
        pvalue_A_vs_B = pvalue(A_count, B_count, library_scaling_factor, pseudo_count)
        pvalue_A_vs_B_list.append(pvalue_A_vs_B)
        pvalue_B_vs_A = pvalue(B_count, A_count, 1 / library_scaling_factor, pseudo_count)
        pvalue_B_vs_A_list.append(pvalue_B_vs_A)

    # Calculate the FDR
    fdr_A_vs_B_list = fdr(pvalue_A_vs_B_list)
    fdr_B_vs_A_list = fdr(pvalue_B_vs_A_list)

    # Output the islands read counts, normalized read counts, fc, pvalue both ways
    scaling_factor = 1000000
    outfile = open(opt.out_file, 'w')
    outline = '#chrom' + "\t" + 'start' + "\t" + 'end' + "\t" + "Readcount_A" + "\t" + 'Normalized_Readcount_A' + "\t" \
              + 'ReadcountB' + "\t" + 'Normalized_Readcount_B' + "\t" + "Fc_A_vs_B" + "\t" + "pvalue_A_vs_B" + "\t" + \
              "FDR_A_vs_B" + "\t" + "Fc_B_vs_A" + "\t" + "pvalue_B_vs_A" + "\t" + "FDR_B_vs_A" + "\n"
    outfile.write(outline)
    ii = 0
    for island in islands_list:
        A_count = island['A_count']
        B_count = island['B_count']
        normalized_A = A_count / float(A_library_size) * scaling_factor
        normalized_B = B_count / float(B_library_size) * scaling_factor
        fc_A_vs_B = ((A_count + pseudo_count) * 1.0 / (B_count + pseudo_count)) / library_scaling_factor
        fc_B_vs_A = ((B_count + pseudo_count) * 1.0 / (A_count + pseudo_count)) * library_scaling_factor
        outline = island['region'].chrom + "\t" + str(island['region'].start) + "\t" + str(island['region'].end) + "\t" + str(
                        A_count) + "\t" + str(normalized_A) + "\t" + str(B_count) + "\t" + str(normalized_B) + "\t" + str(
                        fc_A_vs_B) + "\t" + str(pvalue_A_vs_B_list[ii]) + "\t" + str(fdr_A_vs_B_list[ii]) + "\t" + str(
                        fc_B_vs_A) + "\t" + str(pvalue_B_vs_A_list[ii]) + "\t" + str(fdr_B_vs_A_list[ii]) + "\n"
        outfile.write(outline)
        ii += 1

    # Calculate the correlations using normalized read counts
    A_array = ()
    B_array = ()

    A_array = scipy.array(island_A_readcount_list)
    B_array = scipy.array(island_B_readcount_list)

    # Normalization to reads per million
    A_array = A_array / float(A_library_size) * scaling_factor
    B_array = B_array / float(B_library_size) * scaling_factor
    pearson = scipy.stats.pearsonr(A_array, B_array)
    print "Pearson's correlation is: " + str(pearson[0]) + " with p-value " + str(pearson[1])
    spearman = scipy.stats.spearmanr(A_array, B_array)
    print "Spearman's correlation is: " + str(spearman[0]) + " with p-value " + str(spearman[1])




if __name__ == "__main__":
    main(sys.argv)
