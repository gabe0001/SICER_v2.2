#!/bin/bash
# 
# Authors: Chongzhi Zang, Weiqun Peng, Nick Waring, Ujwal Boddeti
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 2.0  7/24/2017

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

PATHTO=~/SICER2.3
SICER=$PATHTO/SICER

if [ $# -lt 10 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["InputDir"] ["bed file"] ["OutputDir"] ["species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["E-value"] 
    echo ""
    exit 1
fi

#TEMP=$[$3 % $2] 
TEMP=`expr $9 % $6`

if [ $TEMP != 0 ]; then
	echo ""
	echo "Gap_size needs to be multiples of window_size, ie, 0, 2*window_size, etc "
	echo ""
	exit 1
fi

#================================================================================
#Parameters for running SICER without control library and using random background model.

# Input directory
DATADIR=$1

# Input sample bed file
SAMPLEBED=$2

#Output directory
OUTPUTDIR=$3

# Species
SPECIES=$4
SPECIES=${SPECIES:=mm9}

# THRESHOLD is the threshold is for redundancy for reads. THRESHOLD=n
# means that each read has at most n copy after preprocessing.
THRESHOLD=$5
THRESHOLD=${THRESHOLD:=1}

# WINDOW_SIZE is the size of the windows to scan the genome width. 
# One WINDOW_SIZE is the smallest possible island.
WINDOW_SIZE=$6

# FRAGMENT_SIZE is for determination of the amound of shift for center
# of the DNA fragments represented by the reads. FRAGMENT_SIZE=150
# means the shift is 75.
FRAGMENT_SIZE=$7
FRAGMENT_SIZE=${FRAGMENT_SIZE:=150}

# Effective Genome as fraction of the genome size. It depends on read length.
EFFECTIVEGENOME=$8
EFFECTIVEGENOME=${EFFECTIVEGENOME:=0.74}

if [ $(echo "$EFFECTIVEGENOME > 1"|bc) -eq 1 ]; then
	echo ""
	echo " $EFFECTIVEGENOME needs to be between 0 and 1 "
	echo "" 
	exit 1
fi

#GAP_SIZE is in base pairs.
GAP_SIZE=$9

#EVALUE is the number of islands expected in random background. The
#EVALUE is used to determine the score threshold s_T.
EVALUE=${10}



result=`python $PATHTO/src/bam_or_bed.py -f $SAMPLEBED`

echo $result

if [ "$result" == "error" ]; then
echo " "
echo "Error: input file $SAMPLEBED must have either .bam or .bed extension"
exit 1
fi

echo "#############################################"
echo "######           SICER v2.0            ######"
echo "#############################################"

echo "Input library directory: $DATADIR"
echo "ChIP library: $SAMPLEBED"
echo "Output directory: $OUTPUTDIR"
echo "Species: $SPECIES"
echo "Threshold for redundancy allowed for reads: $THRESHOLD"
echo "Window size: $WINDOW_SIZE bps"
echo "Fragment size: $FRAGMENT_SIZE bps. The shift for reads is half of $FRAGMENT_SIZE"
echo "Effective genome size as a fraction of the reference genome of $SPECIES: $EFFECTIVEGENOME"
echo "Gap size: $GAP_SIZE bps"
echo "Evalue for identification of significant islands: $EVALUE"
#================================================================================


if [ "$result" == "bed" ]; then
echo " "
echo " "
echo "python $PATHTO/src/SICER-rb.py -b $SAMPLEBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -e $EVALUE -i $DATADIR -o $OUTPUTDIR -a $PATHTO"
python $PATHTO/src/SICER-rb.py -b $SAMPLEBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -e $EVALUE -i $DATADIR -o $OUTPUTDIR -a $PATHTO

echo "Done!"
fi

if [ "$result" == "bam" ]; then
echo " "
echo " "
echo "python $PATHTO/src/SICER-rb_bam.py -b $SAMPLEBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -e $EVALUE -i $DATADIR -o $OUTPUTDIR -a $PATHTO"
python $PATHTO/src/SICER-rb_bam.py -b $SAMPLEBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -e $EVALUE -i $DATADIR -o $OUTPUTDIR -a $PATHTO

echo "Done!"
fi

