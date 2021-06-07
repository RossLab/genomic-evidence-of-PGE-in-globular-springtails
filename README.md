# Genomic evidence of paternal genome elimination in globular springtails.

This is a supplementary material. It allows, nearly perfect, replication of the study and exact code that was used to generate all the results and figures.

### Springtails collected and sequenced

The sequenced materials resequencing data: PRJEB44694

#### organisation of the input data

`data/reseq/<ID>/...`
`data/reference/<ID>/...`


#### reference genome and X and A assignments



### Coverage estimates

#### k-mer based

#### mapping based

The reads are mapped on the reference using `bowtie2` and sorted using `samtools` as follows

```sh
bowtie2 --very-sensitive-local -p 16 -x $REFERENCE \
        -1 $R1 -2 $R2 \
        --rg-id "$RG_ID" --rg SM:"$SAMPLE" --rg PL:ILLUMINA --rg LB:LIB_"$SAMPLE" \
        | samtools view -h -q -20 \
        | samtools sort -@10 -O bam - > $SAMPLE.rg.sorted.bam
```

**Allacma**

Generating the reference

TODO

Mapping of the reference reads to reference genome

```
scripts/mapping_reference_reads_Afus1.sh
```

**Dicyrtomina**



### Expected coverages of heterozygous loci

Calculation of expectation of allele coverages under different scenarios
    Table 1: table of expected coverages in the two males

### variant calling

Allele coverage distributions


### hypothesis testing
