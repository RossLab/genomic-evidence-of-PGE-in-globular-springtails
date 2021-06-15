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

- `data/raw_reads/{individual}/{accesion}_1.fastq.g`
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
```

**Dicyrtomina**

```
samtools depth data/mapped_reads/Dorn_HP1-3_reads2asm0.bam | scripts/depth2depth_per_contig_median.py > data/mapped_reads/per_scf_cov_medians_Dorn_HP1-3.tsv
samtools depth data/mapped_reads/Dorn_BH2-1_reads2asm0.bam | scripts/depth2depth_per_contig_median.py > data/mapped_reads/per_scf_cov_medians_Dorn_BH2-1.tsv
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

qsub -o logs -e logs -cwd -N SNP_filt -V -pe smp64 1 -b yes 'scripts/subset_SNPs_to_chromosomes.sh data/SNP_calls/freebayes_Ocin2_raw.vcf tables/chr_assignments_Afus1.tsv'
```

```
scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_Afus_filt_sorted data/SNP_calls/freebayes_Afus_filt.vcf tables/chr_assignments_Afus1.tsv

scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_BH2-1_filt_sorted data/SNP_calls/freebayes_BH2-1_filt.vcf tables/chr_assignments_Dorn_BH2-1.tsv

scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_HP1-3_filt_sorted data/SNP_calls/freebayes_HP1-3_filt.vcf tables/chr_assignments_Dorn_HP1-3.tsv

scripts/sort_variants_with_respect_to_chromosomes.py -o data/SNP_calls/freebayes_Ocin2_filt_sorted data/SNP_calls/freebayes_Ocin2_filt.vcf tables/chr_assignments_Ocin1.tsv
```

### hypothesis testing

```
data/SNP_calls/freebayes_all_samples_raw_filt_asn_only.vcf
data/SNP_calls/freebayes_Ocin2_raw_filt_asn_only.vcf
```

### Heterozygous SNP analysis


First, I will extract just ne needed info about the variants in BH3-2

```bash
i=10
for ind in WW5-6 WW5-4 BH3-2 WW3-1 WW2-1 WW1-2 WW2-5 WW5-1 WW5-5 WW1-4 WW2-6 WW5-3; do
  echo $ind $i;
  # grep "^#" data/SNP_calls/freebayes_all_samples_filt_asn_only.vcf | tail -1 | cut -f 1,2,$i
  grep -v "#" data/SNP_calls/freebayes_all_samples_filt_asn_only.vcf | cut -f 1,2,$i | tr -s ":" "\t" | cut -f 1,2,3,4,6,8 | cut -f 1 -d ',' | awk '{ if ($4 != "."){ print $0 } }' > data/SNP_calls/"$ind"_asn_snps.tsv
  i=$(($i+1));
done
```

Now `data/SNP_calls/BH3-2_asn_snps.tsv` contains scf, pos, genotype, total_cov, ref_cov, alt_cov comuns.


```bash
for ind in WW5-6 WW5-4 BH3-2 WW3-1 WW2-1 WW1-2 WW2-5 WW5-1 WW5-5 WW1-4 WW2-6 WW5-3; do
  Rscript scripts/plot_heterozygous_autosomal_allele_coverages.R "$ind" data/SNP_calls/"$ind"_asn_snps.tsv figures/het_autosomal_allele_supports/"$ind".png;
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

# allele_1 <- rnbinom(nrow(male_A_snps), 15, mu = 11.792)
# allele_2 <- rnbinom(nrow(male_A_snps), 15, mu = 11.792)
# total <- allele_1 +  allele_2
#
# minor_cov <- apply(data.frame(cov1 = allele_1, cov2 = allele_2), 1, min)
#
#
# hist(minor_cov, col = 'blue', main = 'Coverages generated by negative binomial and sorted to major/minor', xlab = 'Coverage', breaks = 30)
# hist(total - minor_cov, col = rgb(1, 0.05, 0.1, 0.75), add = T, breaks = 30)

```
