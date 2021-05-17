#!/bin/bash

SAMPLE=Afus1

SCRATCH=/scratch/"$USER"/"$JOB_ID"/mapping
DATADIR=`pwd`
OUTDIR=data/Afus1/mapping
REFERENCE=$DATADIR/data/Afus1/genome/idx

mkdir -p $DATADIR/$OUTDIR
mkdir -p $SCRATCH/$OUTDIR
cd $SCRATCH

R1="data/Afus1/trimmed_reads/Afus1-trimmed-pair1.fastq.gz"
R2="data/Afus1/trimmed_reads/Afus1-trimmed-pair1.fastq.gz"

# handling readgroups from sequencing files
FLOWCELL="HHMYHCCXY"
LANE="008"
RG_ID="$SAMPLE"_"$FLOWCELL"_"$LANE"

bowtie2 --very-sensitive-local -p 16 -x $REFERENCE \
        -1 $R1 -2 $R2 \
        --rg-id "$RG_ID" --rg SM:"$SAMPLE" --rg PL:ILLUMINA --rg LB:LIB_"$SAMPLE" \
        | samtools view -h -q -20 \
        | samtools sort -@10 -O bam - > $SAMPLE.rg.sorted.bam

rsync -av --remove-source-files $SAMPLE.rg.sorted.bam $DATADIR/$OUTDIR

rm BAM_FILES
find /scratch/$USER/$JODID -depth -type d -exec rmdir {} \;
