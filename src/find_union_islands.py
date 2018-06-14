#Ujwal Boddeti - 20 June 2017
#Island Union Program

#Program aims to take two .bed files as inputs and produce an output file which contains
#the union of the islands in both bed files

import os
import sys
import HTSeq
from optparse import OptionParser

def main(argv):
    #establish options for running program
    parser = OptionParser()
    parser.add_option("-a", "--islandfile1", action="store", type="string", dest="file1", metavar="<file>",
                      help="name of first islands file")
    parser.add_option("-b", "--islandfile2", action="store", type="string", dest="file2", metavar="<file>",
                      help="name of second islands file")
    parser.add_option("-o", "--outputfile", action="store", type="string", dest="outfile", metavar="<file>",
                      help="output file name")
    (opt, args) = parser.parse_args(argv)

    print "########## Island Union Program ##########"

    #add file1 islands to tempfile
    os.system('cat %s > %s' % (opt.file1, "tempfile.bed"))

    #add file2 islands to tempfile
    os.system('cat %s >> %s' % (opt.file2, "tempfile.bed"))

    #sort tempfile and store in sortedfile
    os.system('sort -k1,1 -k2,3n %s > %s' % ("tempfile.bed", "sortedfile.bed"))

    #instantiate HTSeq bediterator
    bed_iterator = HTSeq.BED_Reader("sortedfile.bed")

    outfile = open(opt.outfile, 'w')

    total_islands = 0
    tempGI = None
    currentGI = None
    #iterate through GenomicInterval objects
    for read in bed_iterator:
        if tempGI is None:
            currentGI = read.iv
            tempGI = read.iv
        else:
            tempGI = currentGI
            currentGI = read.iv

            #use genomicInterval overlaps method
            if tempGI.overlaps(currentGI):
                currentGI.extend_to_include(tempGI)

            else:
                newLine = str(tempGI.chrom) + "\t" + str(tempGI.start) + "\t" + str(tempGI.end) + "\n"
                outfile.write(newLine)
                total_islands += 1

    #add last entry to union file
    newLine = str(currentGI.chrom) + "\t" + str(currentGI.start) + "\t" + str(currentGI.end) + "\n"
    outfile.write(newLine)
    total_islands += 1

    outfile.close()

    #remove tempfile.bed and sortedfile.bed
    os.system('rm %s %s' % ("tempfile.bed", "sortedfile.bed"))

    print "Total number of islands in islands_union_file: " + str(total_islands)


if __name__ == "__main__":
    main(sys.argv)

#Ujwal Boddeti - 20 June 2017
