#!/bin/bash

READS_DIR=$1
SAMPLE=$2
OUT=$3

SCRATCH=/scratch/"$USER"/"$JOB_ID"/mapping
DATADIR=`pwd`
OUTDIR=`dirname $OUT`
REFERENCE=$DATADIR/data/Afus1/genome/idx

mkdir -p $DATADIR/$OUTDIR
mkdir -p $SCRATCH/$OUTDIR
cd $SCRATCH

for R1 in $DATADIR/$READS_DIR/*-trimmed-pair1.fastq.gz; do
  R2=${R1%1.fastq.gz}2.fastq.gz
  BASE=$(basename $R1 -trimmed-pair1.fastq.gz)

  # handling readgroups from sequencing files
  FLOWCELL=$(echo $BASE | cut -f 3 -d _)
  LANE=$(echo $BASE | cut -f 4 -d _)
  RG_ID="$SAMPLE"_"$FLOWCELL"_"$LANE"

  bowtie2 --very-sensitive-local -p 16 -x $REFERENCE \
          -1 $R1 -2 $R2 \
          --rg-id "$RG_ID" --rg SM:"$SAMPLE" --rg PL:ILLUMINA --rg LB:LIB_"$SAMPLE" \
          | samtools view -h -q -20 \
          | samtools sort -@10 -O bam - > $BASE.bam

  echo $BASE.bam >> BAM_FILES

done

if [ $(cat BAM_FILES | wc -l) -gt 1 ]; then
  samtools merge -b BAM_FILES $SAMPLE.rg.sorted.bam && rm $(cat BAM_FILES)
else
  mv $BASE.bam $SAMPLE.rg.sorted.bam
fi

rsync -av --remove-source-files $SAMPLE.rg.sorted.bam $DATADIR/$OUTDIR

rm BAM_FILES
find /scratch/$USER/$JODID -depth -type d -exec rmdir {} \;
