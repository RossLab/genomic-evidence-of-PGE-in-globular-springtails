#!/bin/bash

R1=$1
R2=$2
OUT=$3

SCRATCH=/scratch/"$USER"/"$JOB_ID"/trim_reads
DATADIR=`pwd`
RAW_READS_DIR=`dirname $R1`
OUTDIR=`dirname $OUT`

mkdir -p $DATADIR/$OUTDIR
mkdir -p $SCRATCH/$OUTDIR
cd $SCRATCH

skewer -z -m pe -n -q 26 -l 21 -t 16 -o $OUT $DATADIR/$R1 $DATADIR/$R2

rsync -av --remove-source-files $SCRATCH/$OUTDIR/* $DATADIR/$OUTDIR

# remove the directory ceated for this job ID to clear up any unwanted files
#rm -rf /scratch/$USER/$JOB_ID

find /scratch/$USER/$JODID -depth -type d -exec rmdir {} \;
