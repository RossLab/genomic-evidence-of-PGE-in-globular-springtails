# Genomic evidence of paternal genome elimination in globular springtails.

This is a supplementary material. It allows, nearly perfect, replication of the study and exact code that was used to generate all the results and figures.

You can also check [ESEB 2021 talk about preliminary results of this work](https://youtu.be/LKBG5AqkUqg?t=4646).

## Springtails collected and sequenced

The sequenced materials resequencing data: PRJEB44694

### organisation of the data

All raw reads will be automatically downloaded from EBI using

```
snakemake download_all_reads
```

Downloading the reference genome

```
mkdir -p data/reference/Afus1/ data/reference/Ocin1/
wget ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/caj/CAJVCH01.fasta.gz -O data/reference/Afus1/genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lji/LJIJ01.fasta.gz -O data/reference/Ocin1/genome.fa.gz
```

The data

**Allacma fusca**

- `data/raw_reads/{individual}/{accesion}_1.fastq.gz`
- `data/reference/Afus1/genome.fa.gz`


## Coverage estimates

This section shows

1. k-mer based estimate of coverage (k-mer coverage) using genomescope
2. mapping of reads and subsequent estimate of mapping coverages

#### k-mer based estimation of two tissue model

This approach is complementary to to the mapping approach, the advantage is that it can be done directly on the k-mer spectra and it's does not require assembling and mapping reads (and therefore avoids all the problems with these processes). The drawback is that each male k-mer spectrum will have co-founded paternal heterozygous sites with X chromosome and maternal autosomal sites, therefore the interpretation of the first k-mer peak is not as straightforward as of the first peak of the mapping coverage.

#### The choice of k

Lower k helps with low coverage datasets in cost of resolution. Usually for genomes with relatively modest coverage and short reads a default k = 21 is used. We also use this k for all our genomes with exception of BH3-2, where we (after attempting of fitting a model with k = 21) reduced the k to 17.

With k = 17, a smaller fraction of genome occurs on unique k-mers, but the k-mer coverage is (R - k + 1 / R) * C, hence lower k is, higher coverage we get. The improvement is only marginal, but we needs only a marginal improvement to make the fit easier.

#### Generating k-mer spectra histograms

Here is an exmaple for one sample (BH3-2) and one values of k (21) - funny enough, this not a k-mer spectrum we used in the end as we later redid it with k = 17, but the principle remains the same.

We start by creating kmer database and then extracting a k-mer histogram

```bash
mkdir data/raw_reads/BH3-2/tmp
ls data/raw_reads/BH3-2/ERR5883998_[1,2].fastq.gz > data/raw_reads/BH3-2/FILES
# kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 10000x
kmc -k21 -t30 -m64 -ci1 -cs100000 @data/raw_reads/BH3-2/FILES data/BH3-2/kmer_db/kmer_counts data/BH3-2/tmp
kmc_tools transform data/BH3-2/kmer_db/kmer_counts histogram data/Ocin1/BH3-2/kmer_k21.hist -cx100000
```

And do this for all the samples.

#### estimating 1n and 2n coverage from k-mer spectra

These models are largely inspired by [GenomsScope source code](https://github.com/schatzlab/genomescope) and GenomeScope is included among the used models. They are adjusted and simplified models that are estimating 1n and 2n peaks separately (i.e. approach similar to the coverage analysis of mapped reads). The simplification of GenomeScope model lays is two aspects - first, we completely ignore anything that occurs in the genome more than once or twice (k-mers in the 1n and 2n peaks), compared to GenomeScope which has an explicit fit of duplications within the genome. Second is that we removed the term that is estimating genome-wise heterozygosity knowing k-mer length. Instead we use term `het`, which simply proportion of 1n and 2n k-mers that are in the 1n peak. That is because we don't actually desire to estmate heterozygosity and the model converges better we simplify the term. The trick to model the two peaks separatelly is in adding of a term `proportion` that models the coverage ratio between first and the second peak, which GenomeScope has fixed to 1:2.

```{R}
# model, x - coverages, y - coverage frequencies, kmerEst, lengthEst, hetEst - starting values for the respective parameters
# note the starting value for proportion is 0.5, meaning that out starting values is the GenomeScope-like estimate
nlsLM_2peak_proportional_peaks <- function(x, y, kmerEst, lengthEst, hetEst = 0.6){
  two_tissue_model <- nlsLM(y ~ ((het       * dnbinom(x, size = (proportion * kmercov2)  / bias, mu = (proportion * kmercov2))) +
                                ((1 - het)  * dnbinom(x, size = kmercov2  / bias, mu = kmercov2))) * length,
          start = list(proportion = 0.5, kmercov2 = 2 * kmerEst, bias = 0.5, length = lengthEst, het = hetEst),
          control = list(minFactor=1e-12, maxiter=40))
  return(two_tissue_model)
}
```

All the fits of the two tissue model are calculated in the [scripts/kmer_estimates_of_two_tissue_model.R](scripts/kmer_estimates_of_two_tissue_model.R) script.

The main reason we formulate the model in this way is that we can calculate asymptotic confidence intervals of all fitted variables, including `proportion`. Now, if the CI of `proportion` spans 0.5, we conclude there was no significant deviation from 1:2 coverage ratio.

To run the model on all the calculated k-mer spectra, run

```
Rscript scripts/kmer_estimates_of_two_tissue_model.R
```

This will generate a bunch of plots and a [table](tables/two_tissue_model_empirical.tsv) with estimates for 2 _A. fusca_ males, 3 _A. fusca_ females and one _O. cincta_ (non-PGE springtial) male (the table is included).  

To test how the method works on samples with various sizes of X chromosomes, levels of heterozygosity, fraction of sperm and various sequencing coverages, we performed a **Power analysis**, check the [documentation](misc/power_analysis_for_estimating_fraction_of_sperm.md) for more details. The summary table with all the results is [here](tables/power_analysis_complete.tsv), to plot the analyis for all the scenarios, run

```
Rscript scripts/plot_power_analysis.R
```

We have also used classical GenomeScope model to validate our approach. The genomescope was fit as follows

```bash
genomescope.R -i data/Ocin1/kmer_db/kmer_k21.hist -o data/Ocin1/genome_profiling -n Ocin1
```

### mapping based coverage estimates

The reads are mapped on the reference using `bowtie2` and sorted using `samtools` as follows

```sh
# indexing reference
# bowtie2-build data/reference/Afus1/genome.fa.gz $REFERENCE

bowtie2 --very-sensitive-local -p 16 -x $REFERENCE \
        -1 $R1 -2 $R2 \
        --rg-id "$RG_ID" --rg SM:"$SAMPLE" --rg PL:ILLUMINA --rg LB:LIB_"$SAMPLE" \
        | samtools view -h -q -20 \
        | samtools sort -@10 -O bam - > $SAMPLE.rg.sorted.bam
```

for all the samples run

```
snakemake map_all
```

Note, when I was mapping reads, I messed up something with the reference sample and run this script `scripts/mapping_reference_reads_Afus1.sh` instead. However, I think I fixed the problem afterwards so this note should be (hopefully) outdated.

**Allacma**

Generate per scaffold coverages

```
samtools depth data/mapped_reads/BH3-2.rg.sorted.rmdup.bam | scripts/depth2cov_per_contig.py > data/mapped_reads/BH3-2_cov_per_scf.tsv
samtools depth data/mapped_reads/Afus1.rg.sorted.rmdup.bam | scripts/depth2cov_per_contig.py > data/mapped_reads/Afus1_cov_per_scf.tsv
```

etc. We then pull together all the mean coverages in a table [tables/Afus_mean_coverage_table.tsv](tables/Afus_mean_coverage_table.tsv) provided in this repo.

**Orchesella**

```
samtools depth data/mapped_reads/Ocin2.rg.sorted.bam | scripts/depth2cov_per_contig.py > data/mapped_reads/Ocin2_cov_per_scf.tsv
```

#### Sexing samples / estimating mapping coverages

We use weighted kernel smoothing to see how many modes there are in the mapping coverage (1 = female, 2 = male). Following script will plot `figures/Allacma_cov_estimates.png` and a bunch of thers too.

```
Rscript scripts/plot_mapping_coverages.R
```

but most importandly, it generates [tables/resequencing_coverage_estimates.tsv](tables/resequencing_coverage_estimates.tsv), which contains estimates of 2n mapping coverage estimates for all the samples, but also 1n estimates for the two male samples (`BH3-2` and `Afus1`).

#### Assigning X chromosomes

To generate `tables/chr_assignments_Afus1.tsv` table, run following R script

```
Rscript scripts/assign-X-linked-scaffolds.R
```

### Expected coverages of heterozygous loci

Calculation of expectation of allele coverages under different scenarios.

#### k-mer based expecations

The k-mer based expecation of maternal and paternal coverage is tricky. We could - in theory - use a formula to convert kmer coverage to genomic coverage, that is however in practice very imprecise. Instead, we use the estimated fraction of sperm, mapping coverage estimate and the formula we used for the two tissue model.

The original formula is `f_h = 1 - (c_A - c_X) / c_X`, where `f_h` is fraction of sperm, and `c_` are respective coverages of X chromosomes and autosomes. With a few modifications we can adjust the equation to be in a form `c_X = c_A / (2 - f_h)` which is also the expectation for the maternal autosomal alleles. Therefore the paternal expecation is `c_A - c_X`.

In the case of BH3-2 `f_h = 0.353` and `c_A = 28.73` (which I took as mean of all homozygous variants on autosomes so we avoid any circularity in the argument). Then `c_m = c_X = 17.44` and therefore `c_p = 11.29`.

#### Coverage estimates

We ploted the modality of the coverage distributions used for sexing of individuals (unimodal - females; bimodal - males; SM Figure 1) and estimated the 1n/2n coverages from mapped reads to the reference using kernel smoothing. Both operation performed in the [plot_mapping_coverages.R](scripts/plot_mapping_coverages.R) script.

As the result, we have [the table](https://github.com/RossLab/genomic-evidence-of-PGE-in-globular-springtails/blob/master/tables/resequencing_coverage_estimates.tsv) with coverage estimates.

### Testing PGE using allele coverages

To test PGE, we called variants and sort them out to chromosomes. The test was based on expectation that the sperm in male _A. fusca_ contains the maternal genome only. If that's the case, there is an expected coverage shift of maternal to paternal alleles which is proportional to the fraction of sperm in the sample (see calculated expectation above).

#### variant calling

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

First, I will extract just the needed info about the variants in BH3-2, Afus1 and Ocin2

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

#### hypothesis testing

Now we have the SNP calls, we need to sort out if one allele have always a higher coverage than the other and that this difference is not simply due to coverage variation.

### Heterozygous SNP analysis

```bash
for ind in WW5-6 WW5-4 BH3-2 WW3-1 WW2-1 WW1-2 WW2-5 WW5-1 WW5-5 WW1-4 WW2-6 WW5-3; do
  Rscript scripts/plot_heterozygous_autosomal_allele_coverages.R "$ind" data/SNP_calls/run0/"$ind"_asn_snps.tsv figures/het_autosomal_allele_supports/"$ind";
done
```

To plot autosomal and X chromosome covera supports (`figures/het_autosomal_allele_supports/autosomal_and_X_variant_coverages_<SAMPLE>.png`) run to get the male plots

```bash
Rscript scripts/plot_heterozygous_autosomal_allele_coverages_BH3-2.R
Rscript scripts/plot_heterozygous_autosomal_allele_coverages_Afus1.R
Rscript scripts/plot_heterozygous_autosomal_allele_coverages_Ocin2.R
```

and for females run

All females

```bash
  for ind in WW5-6 WW5-4 BH3-2 WW3-1 WW2-1 WW1-2 WW2-5 WW5-1 WW5-5 WW1-4 WW2-6 WW5-3; do
   Rscript scripts/plot_heterozygous_autosomal_allele_coverages.R "$ind" data/SNP_calls/run0/"$ind"_asn_snps.tsv figures/het_autosomal_allele_supports/"$ind".png;
  done
```

### data to be deposited:

Heterozygous SNP analysis:
 - `data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_X.tsv`
 - `data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_A.tsv`
 - `data/mapped_reads/Afus1_per_scf_cov_medians.tsv`
 - `data/SNP_calls/freebayes_Afus_filt_sorted_Afus1_A.tsv`
 - `data/SNP_calls/freebayes_Afus_filt_sorted_Afus1_X.tsv`
 - `data/SNP_calls/run0/*_asn_snps.tsv`

Coverage plots:
  - `data/mapped_reads/Afus1_cov_per_scf.tsv`
  - `data/mapped_reads/BH3-2_cov_per_scf.tsv`
  - `data/mapped_reads/Ocin2_cov_per_scf.tsv`

## Supplementary analyses


### Potential causes of the inaccurately high estimate of sperm in Ocin2

Using mapped reads of `Ocin2` to the reference genome (ID), we estimated 1n and 2n coverage to be 53.43x and 101.12x respectively. According to the two tissue model this would imply that `Ocin2` - a species with completely regular meiosis, has ~10% of the body form of cells with the same number of sex chromosomes and

```{R}
asn_tab <- read.table('tables/chr_assignments_Ocin1.tsv', header = T)
ind = 'Ocin2'
ind_tab <- asn_tab[asn_tab[, 'male_coverage'] < filt_quantile & asn_tab[, 'len'] > 20000, c('scf', 'len', 'male_coverage')]

cov_2n <- 101.1228
cov_1n <- 53.4376

png('figures/Ocin2_male_coverage_vs_scf_length.png')
  plot(ind_tab$male_coverage, log10(ind_tab$len), pch = 20, xlab = 'Coverage', ylab = 'log10 Scaffold length', main = 'mapping coverage vs scaffold length in O. cincta male')
  lines(c(53.3, 53.3), c(min(log10(ind_tab$len)), 6), lwd = 2, col = 'red')
  lines(c(101.15, 101.15), c(min(log10(ind_tab$len)), 6), lwd = 2, lty = 2, col = 'red')
  legend('topright', bty = 'n', c('1n coverage est', '2n coverage est'), lty = c(1, 2), lwd = 2, col = 'red')
dev.off()
```

Looking at this plot, it is very clear that at leas some of the long scaffolds are in part missabled or with large scale structural variations. Scaffolds that have a few hundred thousand bases are expected to have the mean mapping coverage very close to the expected value, however, we do see substantial number of data points that are in between of the two distributions and therefore are likely representing a chimera of long stretches of sequences that are in one or two genomic copies in the male genome. These streches might be either large scale deletions or misassembled pieces of the X chromosome.

Every assembly error or structural variation will cause that the mean scaffold mapping covrerage will be somewhere between true 1n and 2n coverage. Hence, the two peaks will always be closer to each other than they would exected to be in perfect data which implies that these errors will bias and inflate the estimate.
