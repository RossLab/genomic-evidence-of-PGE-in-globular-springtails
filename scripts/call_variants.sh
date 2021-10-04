#!/bin/bash

# Afus1
freebayes -f data/reference/Afus1/genome.fa -b data/mapped_reads/Afus1.rg.sorted.rmdup.bam --standard-filters --min-coverage 5 -p 2 > data/SNP_calls/freebayes_Afus1_raw.vcf

# with reference
freebayes -f data/reference/Afus1/genome.fa -b data/mapped_reads/Afus1.rg.sorted.rmdup.bam data/mapped_reads/BH3-2.rg.sorted.rmdup.bam  data/mapped_reads/WW2-1.rg.sorted.rmdup.bam  data/mapped_reads/WW3-1.rg.sorted.rmdup.bam  data/mapped_reads/WW5-4.rg.sorted.rmdup.bam data/mapped_reads/WW1-2.rg.sorted.rmdup.bam  data/mapped_reads/WW2-5.rg.sorted.rmdup.bam  data/mapped_reads/WW5-1.rg.sorted.rmdup.bam  data/mapped_reads/WW5-5.rg.sorted.rmdup.bam data/mapped_reads/WW1-4.rg.sorted.rmdup.bam  data/mapped_reads/WW2-6.rg.sorted.rmdup.bam  data/mapped_reads/WW5-3.rg.sorted.rmdup.bam  data/mapped_reads/WW5-6.rg.sorted.rmdup.bam --populations tables/populations.tsv --hwe-priors-off --standard-filters --min-coverage 5 -p 2 > data/SNP_calls/freebayes_all_samples_raw.vcf

# Ocin2
freebayes -f data/reference/Ocin1/genome.fa -b data/mapped_reads/Ocin2.rg.sorted.bam --standard-filters --min-coverage 5 -p 2 > data/SNP_calls/freebayes_Ocin2_raw.vcf
