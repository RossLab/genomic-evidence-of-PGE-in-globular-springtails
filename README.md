# Genomic evidence of paternal genome elimination in globular springtails.

This is a supplementary material. It allows, nearly perfect, replication of the study and exact code that was used to generate all the results and figures.

You can also check [ESEB talk about preliminary results of this work](https://www.youtube.com/watch?v=LKBG5AqkUqg&t=5363s).

### Springtails collected and sequenced

The sequenced materials resequencing data: PRJEB44694

#### organisation of the data

All data will be automatically downloaded from EBI using

```
./snakemake_clust.sh download_all_reads
```

The data

**Allacma fusca**

- `data/raw_reads/{individual}/{accesion}_1.fastq.gz`
- `data/reference/Afus1/genome.fa.gz`

#### reference genome and X and A assignments



### Coverage estimates

Some numbers

```
# > 95.3122 / 2
# [1] 47.6561
# > 29.6378 / 2
# [1] 14.8189
# > 18.50149 - 29.6378 / 2
# [1] 3.68259
# > 58.09453 - 47.6561
# [1] 10.43843
# > 3.68259 / 29.6378
# [1] 0.1242531
# > 10.43843 / 95.3122
# [1] 0.1095183
```

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

for all the samples run

```
./snakemake_clust.sh map_all
```

**Allacma**

Generating the reference

TODO

Mapping of the reference reads to reference genome

```
scripts/mapping_reference_reads_Afus1.sh

samtools depth data/mapped_reads/BH3-2.rg.sorted.rmdup.bam | scripts/depth2cov_per_contig.py > data/mapped_reads/BH3-2_cov_per_scf.tsv
samtools depth data/mapped_reads/Afus1.rg.sorted.rmdup.bam | scripts/depth2cov_per_contig.py > data/mapped_reads/Afus1_cov_per_scf.tsv
```

**Orchesella**

```
samtools depth data/mapped_reads/Ocin2.rg.sorted.bam | scripts/depth2cov_per_contig.py > data/mapped_reads/Ocin2_cov_per_scf.tsv
```

### Expected coverages of heterozygous loci

Calculation of expectation of allele coverages under different scenarios
    Table 1: table of expected coverages in the two males

### variant calling

We call variants using freebayes.

```
# A. fusca
freebayes -f data/reference/Afus1/genome.fa -b data/mapped_reads/Afus1.rg.sorted.rmdup.bam data/mapped_reads/BH3-2.rg.sorted.rmdup.bam  data/mapped_reads/WW2-1.rg.sorted.rmdup.bam  data/mapped_reads/WW3-1.rg.sorted.rmdup.bam  data/mapped_reads/WW5-4.rg.sorted.rmdup.bam data/mapped_reads/WW1-2.rg.sorted.rmdup.bam  data/mapped_reads/WW2-5.rg.sorted.rmdup.bam  data/mapped_reads/WW5-1.rg.sorted.rmdup.bam  data/mapped_reads/WW5-5.rg.sorted.rmdup.bam data/mapped_reads/WW1-4.rg.sorted.rmdup.bam  data/mapped_reads/WW2-6.rg.sorted.rmdup.bam  data/mapped_reads/WW5-3.rg.sorted.rmdup.bam  data/mapped_reads/WW5-6.rg.sorted.rmdup.bam --populations tables/populations.tsv --hwe-priors-off --standard-filters --min-coverage 5 -p 2 > data/SNP_calls/freebayes_all_samples_raw.vcf

# Ocin2
freebayes -f data/reference/Ocin1/genome.fa -b data/mapped_reads/Ocin2.rg.sorted.bam --standard-filters --min-coverage 5 -p 2 > data/SNP_calls/freebayes_Ocin2_raw.vcf

# Redoing SNP calls for Afus1
freebayes -f data/reference/Afus1/genome.fa -b data/mapped_reads/Afus1.rg.sorted.rmdup.bam --standard-filters --min-coverage 5 -p 2 > data/SNP_calls/freebayes_Afus1_raw.vcf
```

Allele coverage distributions

```bash
# qsub -o logs -e logs -cwd -N SNP_filt -V -pe smp64 1 -b yes '' -> submission on our cluster
scripts/subset_SNPs_to_chromosomes.sh data/SNP_calls/freebayes_all_samples_raw.vcf tables/chr_assignments_Afus1.tsv
scripts/subset_SNPs_to_chromosomes.sh data/SNP_calls/freebayes_Ocin2_raw.vcf tables/chr_assignments_Ocin1.tsv
scripts/subset_SNPs_to_chromosomes.sh data/SNP_calls/freebayes_Afus1_raw.vcf tables/chr_assignments_Afus1.tsv
```

```
data/SNP_calls/freebayes_Afus_filt.vcf
data/SNP_calls/freebayes_Afus1_raw_filt.vcf
data/SNP_calls/freebayes_Ocin2_raw_filt.vcf
```

#### sorting variants to chromosomes

```
python3 scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_Afus_filt_sorted data/SNP_calls/freebayes_Afus_filt.vcf tables/chr_assignments_Afus1.tsv

python3 scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_Ocin2_filt_sorted data/SNP_calls/freebayes_Ocin2_filt.vcf tables/chr_assignments_Ocin1.tsv
```

First, I will extract just ne needed info about the variants in BH3-2, Afus1 and Ocin2

```bash
i=10
for ind in WW5-3 WW5-4 WW5-5 WW3-1 WW2-1 Afus1 BH3-2 WW2-5 WW5-6 WW5-1 WW1-2 WW2-6 WW1-4; do
  echo $ind $i;
  # grep "^#" data/SNP_calls/freebayes_Afus_filt_sorted_A.vcf | tail -1 | cut -f 1,2,$i
  grep -v "#" data/SNP_calls/freebayes_Afus_filt_sorted_A.vcf | cut -f 1,2,$i | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Afus_filt_sorted_"$ind"_A.tsv
  grep -v "#" data/SNP_calls/freebayes_Afus_filt_sorted_X.vcf | cut -f 1,2,$i | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Afus_filt_sorted_"$ind"_X.tsv
  i=$(($i+1));
done

# if to extract just BH3-2
# grep -v "#" data/SNP_calls/freebayes_Afus_filt_sorted_A.vcf | cut -f 1,2,16 | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_A.tsv
# grep -v "#" data/SNP_calls/freebayes_Afus_filt_sorted_X.vcf | cut -f 1,2,16 | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_X.tsv

# Ocin2
grep -v "^#" data/SNP_calls/freebayes_Ocin2_filt_sorted_A.vcf | cut -f 1,2,10 | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Ocin2_filt_sorted_A_Ocin2.tsv
grep -v "^#" data/SNP_calls/freebayes_Ocin2_filt_sorted_X.vcf | cut -f 1,2,10 | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Ocin2_filt_sorted_X_Ocin2.tsv
```

Now we have a bunch of `data/SNP_calls/*_sorted_X.tsv` and `data/SNP_calls/*_sorted_A.tsv` tables that contains scf, pos, genotype, total_cov, ref_cov, alt_cov comuns.

### hypothesis testing

### Heterozygous SNP analysis

```bash
for ind in WW5-6 WW5-4 BH3-2 WW3-1 WW2-1 WW1-2 WW2-5 WW5-1 WW5-5 WW1-4 WW2-6 WW5-3; do
  Rscript scripts/plot_heterozygous_autosomal_allele_coverages.R "$ind" data/SNP_calls/run0/"$ind"_asn_snps.tsv figures/het_autosomal_allele_supports/"$ind".png;
done
```

### plotting notes


Allele cov ratio plots:
  1. BH2-3 (`figures/het_autosomal_allele_supports/autosomal_and_X_variant_coverages_BH3-2.png`) plotted by

  ```bash
  Rscript scripts/plot_heterozygous_autosomal_allele_coverages_BH3-2.R
  ```

  using:
   -  `data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_X.tsv`
   -  `data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_A.tsv`

  2. Afus1

  ```bash
  Rscript scripts/plot_heterozygous_autosomal_allele_coverages_Afus1.R
  ```

  using:
   - `data/mapped_reads/Afus1_per_scf_cov_medians.tsv`
   - `tables/chr_assignments_Afus1.tsv`
   - `data/SNP_calls/freebayes_Afus_filt_sorted_Afus1_A.tsv`
   - `data/SNP_calls/freebayes_Afus_filt_sorted_Afus1_X.tsv`

   3. Ocin2

   ```bash
   Rscript scripts/plot_heterozygous_autosomal_allele_coverages_Ocin2.R
   ```   

   4. All females

   ```bash
   for ind in WW5-6 WW5-4 BH3-2 WW3-1 WW2-1 WW1-2 WW2-5 WW5-1 WW5-5 WW1-4 WW2-6 WW5-3; do
     Rscript scripts/plot_heterozygous_autosomal_allele_coverages.R "$ind" data/SNP_calls/run0/"$ind"_asn_snps.tsv figures/het_autosomal_allele_supports/"$ind".png;
   done
   ```

Coverage plots:
  kmers:
    ```bash
      for dir in data/genome_profiling/*; do
        SAMPLE=$(echo $dir | cut -f 3 -d "/");
        genomescope.R -i $dir/kmer_k21_full.hist -o $dir -n $SAMPLE;
      done
    ```

  mapping:
    Afus1:
     - `data/mapped_reads/Afus1_cov_per_scf.tsv`

    BH3-2:
     - `data/mapped_reads/BH3-2_cov_per_scf.tsv`

    Ocin2:
     - `data/mapped_reads/Ocin2_cov_per_scf.tsv`


    ```
      Rscript scripts/plot_mapping_coverages.R
    ```

    plots them all (`figures/mapping_coverages/*`).
