## Estimating fractions of sperm using k-mer spectra

This approach is complementary to to the mapping approach, the advantage is that it can be done directly on the k-mer spectra and it's does not require assembling and mapping reads (and therefore avoids all the problems with these processes). The downsite is that each male k-mer spectrum will have co-founded paternal heterozygous sites with X chromosome and maternal autosomal sites, therefore the interpretation of the first k-mer peak is not as streightforward as of the first peak of the mapping coverage.

###re-doing k-mer spectra of BH3-2 with lower k

Lower k helps with low coverage datasets in cost of resolution. With k = 17, a smaller fraction of genome occurs on unique k-mers, but the k-mer coverage is (R - k + 1 / R) * C, hence lower k is, higher coverage we get. The improvement is only marginal, but we needs only a marginal improvement to make the fit easier.

```bash
conda activate ntcard
# gospring
qsub -o logs -e logs -cwd -N ntcard -V -pe smp64 16 -b yes 'time ntcard -t 16 -k 17 -c 10000 -o data/resequencing/genome_profiling/BH3-2/BH3-2_k17.hist @data/resequencing/genome_profiling/BH3-2/FILES'
```

### k-mer spectra to 1n and 2n estimates

These models are largely insprired by [GenomsScope source code](https://github.com/schatzlab/genomescope) and GenomeScope is included among the used models. They are adjusted and simplified models that are estimating 1n and 2n peaks indipendently (i.e. approach similar to the coverage analysis of mapped reads). The simplification of GenomeScope model lays is two aspects - first, we completelly ignore anything that occurs in the genome more than once or twice (kmers in the 1n and 2n peaks), while GenomeScope has an explicit fit of duplications within the genome. Second is that we removed the term that is estimating genome-wise heterozygosity knowing k-mer length. Instead we use term `het`, which simply proportion of 1n and 2n k-mers that are in the 1n peak. That is because we don't actually desire to estmate hererozygosity and the model converges better we simplify the term. Hence the model to fit the two peaks indipendently is

```{R}
# model, x - coverages, y - coverage frequencies, kmerEst, lengthEst, hetEst - starting values for the respective parameters
# note that kmerEst (and 2*) is a starting value for both kmercov and kmercov2 - which are the two indipendent estimates of positions of the coverage peaks, bias is the overdispersal term
nlsLM_2peak_unconditional_peaks <- function(x, y, kmerEst, lengthEst, hetEst = 0.6){
  nlsLM(y ~ ((het       * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
            ((1 - het)  * dnbinom(x, size = kmercov2  / bias, mu = kmercov2))) * length,
        start = list(kmercov = kmerEst, kmercov2 = (2 * kmerEst), bias = 0.5, length = lengthEst, het = hetEst),
        control = list(minFactor=1e-12, maxiter=40))
}
```

### female samples

```{R}
require('minpack.lm')

females <- c('WW5-3', 'WW5-5', 'WW2-6')

WW5_3_k21_file = './data/genome_profiling/WW5-3/kmer_k21_full.hist'
WW5_5_k21_file = './data/genome_profiling/WW5-5/kmer_k21_full.hist'
WW2_6_k21_file = './data/genome_profiling/WW2-6/kmer_k21_full.hist'

source('scripts/genome_models.R') # load genome models (see the script for their details)

### female
feAfus1_kmer_spectrum <- read.table(WW5_3_k21_file, col.names = c('coverage', 'frequency'))
feBH32_kmer_spectrum <- read.table(WW5_5_k21_file, col.names = c('coverage', 'frequency'))
female3_kmer_spectrum <- read.table(WW2_6_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:120
plot(frequency ~ coverage, data = feBH32_kmer_spectrum[cov_range, ])

#
female1_model <- nls_4peak(feAfus1_kmer_spectrum[cov_range, 'coverage'], feAfus1_kmer_spectrum[cov_range, 'frequency'], 21, 15, 300e6, 40)
female2_model <- nls_4peak(feBH32_kmer_spectrum[cov_range, 'coverage'], feBH32_kmer_spectrum[cov_range, 'frequency'], 21, 20, 300e6, 40)
female3_model <- nls_4peak(female3_kmer_spectrum[cov_range, 'coverage'], female3_kmer_spectrum[cov_range, 'frequency'], 21, 15, 300e6, 40)

female_parameters <- as.data.frame(rbind(coef(female1_model), coef(female2_model), coef(female3_model)))
female_parameters$ind <- females
rownames(female_parameters) <- females

write.table(female_parameters, 'tables/kmer_spectra_modeling_female_parameters.tsv', quote = F, col.names = T, row.names = F)
# these are practically genomescope estimates, we can also fit there some of the models that fit indipendently 1n and 2n coverages

# FITTING FEMALE PGE MODELS
# # here, using regular GenomeScope estimates of length and heterozygosity as fixed parameters to reduce complexity of the models
female1_PGE_model <- nlsLM_2peak_unconditional_peaks(feAfus1_kmer_spectrum[cov_range, 'coverage'], feAfus1_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)

female2_PGE_model <- nlsLM_2peak_unconditional_peaks(feBH32_kmer_spectrum[cov_range, 'coverage'], feBH32_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)

female3_PGE_model <- nlsLM_2peak_unconditional_peaks(female3_kmer_spectrum[cov_range, 'coverage'], female3_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)
```

### Ocin2 sample

```{R}
Ocin2_k21_file = './data/genome_profiling/Ocin2/kmer_k21_full.hist'
Ocin2_kmer_spectrum <- read.table(Ocin2_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:150

source('scripts/genome_models.R') # load genome models (see the script for their details)

Ocin2_Genomescope <- nls_4peak(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 21, 50, 280e6, 40)

Ocin2_PGE_model <- nlsLM_2peak_unconditional_peaks(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 50, 300e6, 0.3)
```

### Afus1

Load k-mer spectra of the two males

```{R}
# REFERENCE MALE 1
Afus1_k21_file = './data/genome_profiling/Afus1/kmer_k21_full.hist'
Afus1_kmer_spectrum <- read.table(Afus1_k21_file, col.names = c('coverage', 'frequency'))
Afus1_range <- 20:250

plot(Afus1_kmer_spectrum[Afus1_range, 'frequency'] ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'])

k = 21
eyeballed_Afus1_coverage = 50

### original - works
# source('scripts/genome_models.R')
Afus1_genomescope <- nls_4peak(Afus1_kmer_spectrum[Afus1_range, 'coverage'],
                         Afus1_kmer_spectrum[Afus1_range, 'frequency'],
                         k, eyeballed_Afus1_coverage, 300e6)
lines(predict(Afus1_genomescope, response = T) ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'])

Afus1_PGE_model <- nlsLM_2peak_unconditional_peaks(Afus1_kmer_spectrum[Afus1_range, 'coverage'],  Afus1_kmer_spectrum[Afus1_range, 'frequency'], 55, 280e6, 0.5)
lines(predict(Afus1_PGE_model, response = T) ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'], col = 'red', lty = 2, pwd = 2)
```

### BH3-2

```{R}
# MALE 2
# BH3_male_k21_file = './data/genome_profiling/BH3-2/kmer_k21_full.hist'
BH3_male_k17_file = './data/genome_profiling/BH3-2/BH3-2_k17_reduced.hist'

BH32_kmer_spectrum <- read.table(BH3_male_k17_file, col.names = c('coverage', 'frequency'))
BH32_range <- 5:100

plot(BH32_kmer_spectrum[BH32_range, 'frequency'] ~ BH32_kmer_spectrum[BH32_range, 'coverage'])

k = 17
eyeballed_coverage_BH32 = 15

### original - works
# source('scripts/genome_models.R')
BH32_genomescope <- nls_4peak(BH32_kmer_spectrum[BH32_range, 'coverage'],
                              BH32_kmer_spectrum[BH32_range, 'frequency'],
                              k,
                              eyeballed_coverage_BH32,
                              300e6)

BH32_PGE_model <- nlsLM_2peak_unconditional_peaks(BH32_kmer_spectrum[BH32_range, 'coverage'],  BH32_kmer_spectrum[BH32_range, 'frequency'], eyeballed_coverage_BH32, 280e6, 0.5)
lines(predict(BH32_PGE_model, response = T) ~ BH32_kmer_spectrum[BH32_range, 'coverage'], col = 'red', lty = 2, pwd = 2)

# nlstools::confint2(BH32_PGE_model, level = 0.95, method = "asymptotic")
```

### Summary of all the models

```{R}
inds <- c('Afus1', 'BH3-2', 'Ocin2', females)
GenomeScope_models <- list(Afus1_genomescope, BH32_genomescope, Ocin2_Genomescope, female1_model, female2_model, female3_model)
PGE_models <- list(Afus1_PGE_model, BH32_PGE_model, Ocin2_PGE_model, female1_PGE_model, female2_PGE_model, female3_PGE_model)

PGE_1n_coverages <- sapply(PGE_models, function(m){ coef(m)['kmercov'] })
PGE_2n_coverages <- sapply(PGE_models, function(m){ coef(m)['kmercov2'] })

data.frame(ind = inds, 'cov_1n' = PGE_1n_coverages, 'cov_2n' = PGE_2n_coverages)

source('scripts/genome_models.R')
plot_PGE_model(PGE_models[[1]])


```
