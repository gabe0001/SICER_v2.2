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


def main(argv):

    parser = OptionParser()
    parser.add_option("-f", "--file", action="store", type="string", dest="file_name", metavar="<file>",
                      help="name of file")
    parser.add_option("-g", "--file2", action="store", type="string", dest="file_name2", metavar="<file>",
                      help="name of file2")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 1:
        parser.print_help()
        sys.exit(1)

    extension = opt.file_name[-3:]
    extension = str.lower(extension)

    extension2 = opt.file_name2[-3:]
    extension2 = str.lower(extension2)

    if extension == "bed" and extension2 == "bed":
        print "bed"
    elif extension == "bam" and extension2 == "bam":
        print "bam"
    else:
        print "error"


if __name__ == "__main__":
    main(sys.argv)

