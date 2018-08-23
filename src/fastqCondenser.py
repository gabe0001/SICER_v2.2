#Created by: Ujwal Boddeti on 1 August 2017
#Condensed FASTQ Generator

#Program takes a FASTQ file as an input and produces an output FASTQ file with the first specified number of reads
#and some random number of duplicate reads
#The program also specifies how many total reads should be retained for a specified redundancy threshold after reads that
#don't pass the redundncy threshold are removed

import os
import sys
import HTSeq
from random import randint
from random import seed
from optparse import OptionParser

class Read:
    #attribute references
    key = 0
    line1 = ""
    line2 = ""
    line3 = ""
    line4 = ""

    def __init__(self, arg1, arg2, arg3, arg4, arg5):
        self.key = arg5
        self.line1 = str(arg1)
        self.line2 = str(arg2)
        self.line3 = str(arg3)
        self.line4 = str(arg4)

    def __str__(self):
        return "Read number " + str(self.key)

def printRead(outfile, read):

    outfile.write(read.line1)
    outfile.write(read.line2)
    outfile.write(read.line3)
    outfile.write(read.line4)

def numRetainedReads(reads, dupList, redundancythreshold):
    x = reads
    x = x - len(dupList)
    for y in dupList:
        if y >= redundancythreshold:
            x+=redundancythreshold
        else:
            x+=y
    return x

def main(argv):
    #establish options for running program
    parser = OptionParser()
    parser.add_option("-f", "--fastq", action="store", type="string", dest="file1", metavar="<file>",
                      help="name of bed file")
    parser.add_option("-n", "--reads", action="store", type="int", dest="numReads", metavar="<int>",
                      help="number of desired reads in output file")
    parser.add_option("-c", "--duplicationchance", action="store", type="int", dest="dupchance", metavar="<int>",
                      help="chance a read is duplicated in percent without percent symbol (default is 20)")
    parser.add_option("-d", "--duplicationrange", action="store", type="int", dest="numdups", metavar="<int>",
                      help="range of number of times a read maybe duplicated (default is 10)")
    parser.add_option("-r", "--redundancythreshold", action="store", type="int", dest="redundancythreshold", metavar="<int>",
                      help="how many duplicates want to be retained")
    parser.add_option("-s", "--randomseed", action="store", type="int", dest="randomseed", metavar="<int>",
                      help="random number seed")
    (opt, args) = parser.parse_args(argv)

    infile = open(opt.file1, 'r')
    outfile = open(opt.file1[:-6]+"_minimized"+"-"+str(opt.numReads)+".fastq", 'w')
    outfile1 = open("redundancyTrashCheck.fastq", 'w')
    #outfile2 = open("dupCount.txt", 'w')
    #outfile3 = open("100dupedReads.fastq", 'w')

    numReads = opt.numReads

    count = 0
    reads = 0
    totalReads = 0

    dupchance = 0
    numdups = 0

    dupList = []
    redundancythreshold = opt.redundancythreshold

    key = 0

    #default chance a read is duplicated is 20%
    if opt.dupchance is None:
        dupchance = 5
    else:
        dupchance = 100/opt.dupchance

    #default maximum number of times a read can be duplicated
    if opt.numdups is None:
        numdups = 10
    else:
        numdups = opt.numdups

    print "Number of unique reads in generated FASTQ file: " + str(numReads)

    while reads < numReads:
        #reset random parameters
        seed(opt.randomseed*reads)
        dupChance = randint(1,dupchance)
        numDups = randint(1, numdups)

        read = Read(infile.next(), infile.next(), infile.next(), infile.next(), key)

        #read is written to outfile
        printRead(outfile, read)
        totalReads+=1

        if dupChance == 1:
            dupList.append(numDups+1)
            #printRead(outfile1, read)
            #outfile2.write(str(numDups+1)+"\n")
            while numDups > 0:
                printRead(outfile, read)
                totalReads+=1
                numDups-=1

        key+=1
        count+=4
        reads+=1

    print "Total number of reads: " + str(totalReads) + " (total number of lines is " + str(totalReads*4) + ")"
    print "Number of reads mulitplied is: " + str(len(dupList))
    print "Specified redundancy threshold: " + str(redundancythreshold) + " (how many multiplied reads to keep)"
    print "Number of desired retained reads: " + str(numRetainedReads(reads, dupList, redundancythreshold)) + " (not accounting for unmapped reads)"

    infile.close()
    outfile.close()

if __name__ == "__main__":
    main(sys.argv)

#Ujwal Boddeti - 1 August 2017
