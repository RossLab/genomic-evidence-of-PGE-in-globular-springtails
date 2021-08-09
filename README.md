# Genomic evidence of paternal genome elimination in globular springtails.

This is a supplementary material. It allows, nearly perfect, replication of the study and exact code that was used to generate all the results and figures.

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

samtools depth data/mapped_reads/BH3-2.rg.sorted.rmdup.bam | scripts/depth2depth_per_contig_median.py > data/mapped_reads/per_scf_cov_medians_BH3-2.tsv
samtools depth data/mapped_reads/BH3-2.rg.sorted.rmdup.bam | scripts/depth2depth_per_contig_median.py > data/mapped_reads/per_scf_cov_medians_Afus1.tsv
# scripts/bam2per_scf_medians.sh data/mapped_reads/Afus1.rg.sorted.rmdup.bam data/mapped_reads/Afus1_per_scf_cov_medians.tsv
```

**Orchesella**

```
samtools depth data/mapped_reads/Ocin2.rg.sorted.bam | scripts/depth2depth_per_contig_median.py > data/mapped_reads/per_scf_cov_medians_Ocin2.tsv
```

### Expected coverages of heterozygous loci

Calculation of expectation of allele coverages under different scenarios
    Table 1: table of expected coverages in the two males

### variant calling

Allele coverage distributions

```bash
qsub -o logs -e logs -cwd -N SNP_filt -V -pe smp64 1 -b yes 'scripts/subset_SNPs_to_chromosomes.sh data/SNP_calls/freebayes_all_samples_raw.vcf tables/chr_assignments_Afus1.tsv'

qsub -o logs -e logs -cwd -N SNP_filt -V -pe smp64 1 -b yes 'scripts/subset_SNPs_to_chromosomes.sh data/SNP_calls/freebayes_Ocin2_raw.vcf tables/chr_assignments_Ocin1.tsv'

qsub -o logs -e logs -cwd -N SNP_filt -V -pe smp64 1 -b yes 'scripts/subset_SNPs_to_chromosomes.sh data/SNP_calls/freebayes_Afus_Afus1_raw.vcf tables/chr_assignments_Afus1.tsv'
```

```
scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_Afus_filt_sorted data/SNP_calls/freebayes_Afus_filt.vcf tables/chr_assignments_Afus1.tsv

scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_Ocin2_filt_sorted data/SNP_calls/freebayes_Ocin2_filt.vcf tables/chr_assignments_Ocin1.tsv

scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_Afus1_filt_sorted_by_median data/SNP_calls/freebayes_Afus_Afus1_raw_filt.vcf tables/chr_assignments_Afus1_medians.tsv
```

### hypothesis testing

```
data/SNP_calls/freebayes_all_samples_raw_filt_asn_only.vcf
data/SNP_calls/freebayes_Ocin2_raw_filt_asn_only.vcf
data/SNP_calls/freebayes_Afus_Afus1_raw_filt.vcf
```

### Heterozygous SNP analysis


First, I will extract just ne needed info about the variants in BH3-2

```bash
i=10
for ind in WW5-3 WW5-4 WW5-5 WW3-1 WW2-1 Afus1 BH3-2 WW2-5 WW5-6 WW5-1 WW1-2 WW2-6 WW1-4; do
  echo $ind $i;
  # grep "^#" data/SNP_calls/freebayes_Afus_filt_sorted_A.vcf | tail -1 | cut -f 1,2,$i
  grep -v "#" data/SNP_calls/freebayes_Afus_filt_sorted_A.vcf | cut -f 1,2,$i | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Afus_filt_sorted_"$ind"_A.tsv
  grep -v "#" data/SNP_calls/freebayes_Afus_filt_sorted_X.vcf | cut -f 1,2,$i | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Afus_filt_sorted_"$ind"_X.tsv
  i=$(($i+1));
done

grep -v "^#" data/SNP_calls/freebayes_Ocin2_filt_sorted_A.vcf | cut -f 1,2,10 | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Ocin2_filt_sorted_A_Ocin2.tsv
grep -v "^#" data/SNP_calls/freebayes_Ocin2_filt_sorted_X.vcf | cut -f 1,2,10 | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Ocin2_filt_sorted_X_Ocin2.tsv

grep -v "^#" data/SNP_calls/freebayes_Afus_Afus1_raw_filt.vcf | cut -f 1,2,10 | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($3 != "."){ print $0 } }' > data/SNP_calls/freebayes_Afus_Afus1_filt_reduced.tsv
```

Now `data/SNP_calls/BH3-2_asn_snps.tsv` contains scf, pos, genotype, total_cov, ref_cov, alt_cov comuns.


```bash
for ind in WW5-6 WW5-4 BH3-2 WW3-1 WW2-1 WW1-2 WW2-5 WW5-1 WW5-5 WW1-4 WW2-6 WW5-3; do
  Rscript scripts/plot_heterozygous_autosomal_allele_coverages.R "$ind" data/SNP_calls/run0/"$ind"_asn_snps.tsv figures/het_autosomal_allele_supports/"$ind".png;
done
```

```R
male_snp_tab <- read.table("data/SNP_calls/BH3-2_asn_snps.tsv")
colnames(male_snp_tab) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')

asn_tab <- read.table('tables/chr_assignments_Afus1.tsv', header = T)
strlen <- max(nchar(asn_tab$scf))

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(asn_tab) <- names2tokens(asn_tab$scf)
male_snp_tab$chr <- asn_tab[names2tokens(male_snp_tab$scf), 'chr']

male_snp_tab$cov_minor <- apply(male_snp_tab[, c('ref_cov', 'alt_cov')], 1, min)
male_snp_tab$major_cov <- male_snp_tab$total_cov - male_snp_tab$cov_minor

male_A_snps <- male_snp_tab[male_snp_tab$chr == 'A' & male_snp_tab$genotype == '0/1', ]
male_X_snps <- male_snp_tab[male_snp_tab$chr == 'X' & male_snp_tab$genotype == '1/1', ]

pal <- c('black', rgb(0.8, 0.05, 0.1, 0.55), rgb(0.02, 0.45, 0.65, 0.55))

hist(male_X_snps$total_cov, breaks = 60, main = 'BH2-3', xlab = 'Coverage', freq = F, ylim = c(0, 0.1), col = pal[1])
hist(male_A_snps$major_cov, col = pal[2], add = T, breaks = 120, freq = F)
hist(male_A_snps$cov_minor, col = pal[3], add = T, freq = F, breaks = 60)
legend('topright', c('X 1/1', 'A 0/1 major (maternal)', 'A 0/1 minor (paternal)'), pch = 20, cex = 1.5, bty = 'n', col = pal)
```

Let's try redo `Afus1`

```R
male_snp_tab <- read.table('data/SNP_calls/freebayes_Afus_Afus1_filt_reduced.tsv')
colnames(male_snp_tab) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')

asn_tab <- read.table('tables/chr_assignments_Afus1.tsv', header = T)
strlen <- max(nchar(male_snp_tab$scf))

cov_tab <- read.table("../springtailome/data/resequencing/mapping/complete_coverage_table_medians.tsv", header = T, check.names = F)

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(asn_tab) <- names2tokens(asn_tab$scf)

male_snp_tab$chr <- asn_tab[names2tokens(male_snp_tab$scf), 'chr']

male_snp_tab$cov_minor <- apply(male_snp_tab[, c('ref_cov', 'alt_cov')], 1, min)
male_snp_tab$major_cov <- male_snp_tab$total_cov - male_snp_tab$cov_minor

cov_tab$chr <- asn_tab[names2tokens(cov_tab$scf), 'chr']
male_X <- cov_tab[which(cov_tab$chr == 'X' & cov_tab$chr_medians == 'X'), c('scf', 'len', 'Afus1')]

# male_snp_tab$chr_medians == 'A' &
male_A_snps <- male_snp_tab[which(male_snp_tab$chr == 'A'), ]
male_A_snps <- male_A_snps[male_A_snps$total_cov < 150, ]

pal <- c('black', rgb(0.8, 0.05, 0.1, 0.55), rgb(0.02, 0.45, 0.65, 0.55))

# hist(male_X$Afus1, breaks = 60, main = 'Afus1', xlab = 'Coverage', freq = F, ylim = c(0, 0.1), col = pal[1], xlim = c(0, 140))
hist(male_A_snps$major_cov, col = pal[2], main = 'Afus1', ylim = c(0, 0.06), xlim = c(0, 120), xlab = 'Coverage', breaks = 120, freq = F)
hist(male_A_snps$cov_minor, col = pal[3], add = T, freq = F, breaks = 60)
lines(c(57.45, 57.45), c(0, 100), lwd = 2, lty = 2)

legend('topright', c('X coverage expectation', 'A 0/1 major (maternal)', 'A 0/1 minor (paternal)'), pch = c(NA, 20, 20), lty = c(2, NA, NA), cex = 1.3, bty = 'n', col = pal)
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