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

if [ $# -lt 11 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"]   ["gap size (bp)"] ["FDR"]
    echo ""
    exit 1
fi


TEMP=`expr ${10} % $7`

if [ $TEMP != 0 ]; then
	echo ""
	echo "Gap_size needs to be multiples of window_size, ie, 0, 2*window_size, etc "
	echo ""
	exit 1
fi

#================================================================================
#############################################
# ######  SET UP DIRECTORY STRUCTURE ###### #
#############################################

# The path to the data files.
DATADIR=$1

# Output directory
OUTPUTDIR=$4
#================================================================================

# Input sample bed file
SAMPLEBED=$2

# Input control bed file
CONTROLBED=$3

# Species
SPECIES=$5
SPECIES=${SPECIES:=mm9}

# THRESHOLD is the threshold is for redundancy allowed for reads. THRESHOLD=n
# means that each read has at most n copy after preprocessing.
THRESHOLD=$6
THRESHOLD=${THRESHOLD:=1}

# WINDOW_SIZE is the size of the windows to scan the genome width.
# One WINDOW_SIZE is the smallest possible island.
WINDOW_SIZE=$7

# FRAGMENT_SIZE is for determination of the amound of shift for center
# of the DNA fragments represented by the reads. FRAGMENT_SIZE=150
# means the shift is 75.
FRAGMENT_SIZE=$8
FRAGMENT_SIZE=${FRAGMENT_SIZE:=150}

# Effective Genome as fraction of the genome size. It depends on read length.
EFFECTIVEGENOME=$9
EFFECTIVEGENOME=${EFFECTIVEGENOME:=0.74}

#EVALUE is the number of islands expected in random background. The E value is used for identification of candidate islands that exhibit clustering.
EVALUE=1000

#GAP_SIZE is in base pairs.
GAP_SIZE=${10}

#False discovery rate controlling significance
FDR=${11}

result=`python $PATHTO/src/bam_or_bed_2.py -f $SAMPLEBED -g $CONTROLBED`

if [ "$result" == "error" ]; then
echo " "
echo "Error: input files $SAMPLEBED and $CONTROLBED must have either .bam or .bed extension (and must be the same one)"
exit 1
fi

echo "#############################################"
echo "######           SICER v2.0            ######"
echo "#############################################"

echo "Input library directory: $DATADIR"
echo "ChIP library: $SAMPLEBED"
echo "Control library: $CONTROLBED"
echo "Output directory: $OUTPUTDIR"
echo "Species: $SPECIES"
echo "Threshold for redundancy allowed for chip reads: $THRESHOLD"
echo "Threshold for redundancy allowed for control reads: $THRESHOLD"
echo "Window size: $WINDOW_SIZE bps"
echo "Fragment size: $FRAGMENT_SIZE bps. The shift for reads is half of $FRAGMENT_SIZE"
echo "Effective genome size as a fraction of the reference genome of $SPECIES: $EFFECTIVEGENOME"
echo "Gap size: $GAP_SIZE bps"
echo "Evalue for identification of candidate islands that exhibit clustering: $EVALUE"
echo "False discovery rate controlling significance: $FDR"


if [ "$result" == "bed" ]; then
echo " "
echo " "
echo "python $PATHTO/src/SICER.py -b $SAMPLEBED -c $CONTROLBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -d $FDR -i $DATADIR -o $OUTPUTDIR -a $PATHTO"
python $PATHTO/src/SICER.py -b $SAMPLEBED -c $CONTROLBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -d $FDR -i $DATADIR -o $OUTPUTDIR -a $PATHTO

echo "Done!"
fi

if [ "$result" == "bam" ]; then
echo " "
echo " "
echo "python $PATHTO/src/SICER_bam.py -b $SAMPLEBED -c $CONTROLBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -d $FDR -i $DATADIR -o $OUTPUTDIR -a $PATHTO"
python $PATHTO/src/SICER_bam.py -b $SAMPLEBED -c $CONTROLBED -g $SPECIES -r $THRESHOLD -w $WINDOW_SIZE -f $FRAGMENT_SIZE -p $EFFECTIVEGENOME -s $GAP_SIZE -d $FDR -i $DATADIR -o $OUTPUTDIR -a $PATHTO

echo "Done!"
fi


