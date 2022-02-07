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

These models are largely insprired by [GenomsScope source code](https://github.com/schatzlab/genomescope) and GenomeScope is included among the used models. They are adjusted and simplified models that are either - estimating 1n and 2n peaks indipendently (approach similar to mapping) or fitting 1n coverage of the sperm and soma indipendently (sperm will act as "constant shift" factor).

```

```

### female samples

```{R}
require('minpack.lm')

WW5_3_k21_file = './data/genome_profiling/WW5-3/kmer_k21_full.hist'
WW5_5_k21_file = './data/genome_profiling/WW5-5/kmer_k21_full.hist'
WW2_6_k21_file = './data/genome_profiling/WW2-6/kmer_k21_full.hist'

source('scripts/genome_models.R') # load genome models (see the script for their details)

### female
female1_kmer_spectrum <- read.table(WW5_3_k21_file, col.names = c('coverage', 'frequency'))
female2_kmer_spectrum <- read.table(WW5_5_k21_file, col.names = c('coverage', 'frequency'))
female3_kmer_spectrum <- read.table(WW2_6_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:120
plot(frequency ~ coverage, data = female2_kmer_spectrum[cov_range, ])

#
female1_model <- nls_4peak(female1_kmer_spectrum[cov_range, 'coverage'], female1_kmer_spectrum[cov_range, 'frequency'], 21, 15, 300e6, 40)
female2_model <- nls_4peak(female2_kmer_spectrum[cov_range, 'coverage'], female2_kmer_spectrum[cov_range, 'frequency'], 21, 20, 300e6, 40)
female3_model <- nls_4peak(female3_kmer_spectrum[cov_range, 'coverage'], female3_kmer_spectrum[cov_range, 'frequency'], 21, 15, 300e6, 40)

female_parameters <- as.data.frame(rbind(coef(female1_model), coef(female2_model), coef(female3_model)))
female_parameters$ind <- c('WW5-3', 'WW5-5', 'WW2-6')
rownames(female_parameters) <- c('WW5-3', 'WW5-5', 'WW2-6')

write.table(female_parameters, 'tables/kmer_spectra_modeling_female_parameters.tsv', quote = F, col.names = T, row.names = F)
# these are practically genomescope estimates, we can also fit there some of the models that fit indipendently 1n and 2n coverages

# FITTING FEMALE PGE MODELS
# # here, using regular GenomeScope estimates of length and heterozygosity as fixed parameters to reduce complexity of the models
female1_PGE_model <- nlsLM_2peak_unconditional_peaks(female1_kmer_spectrum[cov_range, 'coverage'], female1_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)

female2_PGE_model <- nlsLM_2peak_unconditional_peaks(female2_kmer_spectrum[cov_range, 'coverage'], female2_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)

female3_PGE_model <- nlsLM_2peak_unconditional_peaks(female3_kmer_spectrum[cov_range, 'coverage'], female3_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)
```

### Ocin2 sample

```{R}
Ocin2_k21_file = './data/genome_profiling/Ocin2/kmer_k21_full.hist'
Ocin2_kmer_spectrum <- read.table(Ocin2_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:150

source('scripts/genome_models.R') # load genome models (see the script for their details)

Ocin2_Genomescope <- nls_4peak(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 21, 50, 280e6, 40)

Ocin2_PGE_model <- nlsLM_2peak_unconditional_peaks(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 50, 280e6, 0.3)
```

### Afus1

Load k-mer spectra of the two males

```{R}
# REFERENCE MALE 1
Afus1_k21_file = './data/genome_profiling/Afus1/kmer_k21_full.hist'
male1_kmer_spectrum <- read.table(Afus1_k21_file, col.names = c('coverage', 'frequency'))
m1_range <- 20:250

plot(male1_kmer_spectrum[m1_range, 'frequency'] ~ male1_kmer_spectrum[m1_range, 'coverage'])

k = 21
eyeballed_coverage_m1 = 50

### original - works
# source('scripts/genome_models.R')
genomescope_m1 <- nls_4peak(male1_kmer_spectrum[m1_range, 'coverage'],
                         male1_kmer_spectrum[m1_range, 'frequency'],
                         k, eyeballed_coverage_m1, 300e6)
lines(predict(genomescope, response = T) ~ male1_kmer_spectrum[m1_range, 'coverage'])

biasEst <- coef(genomescope)['bias']
kmerCoverageEst<- coef(genomescope)['kmercov']
# lengths and genome size better derived from female hist
lengthEst <- mean(female_parameters[, 'length'])
rEst <- mean(female_parameters[, 'r'])

two_tissue_4peak_model <- nls_4peak_two_tissue_constant_shift(male1_kmer_spectrum[m1_range, 'coverage'],
                                                              male1_kmer_spectrum[m1_range, 'frequency'],
                                                              21, 50, lengthEst, 0)
lines(predict(two_tissue_4peak_model, response = T) ~ male1_kmer_spectrum[m1_range, 'coverage'], lwd = 2, col = 'yellow')

two_tissue_model <- nls_2peak_two_tissue_constant_shift(male1_kmer_spectrum[m1_range, 'coverage'],
                                                        male1_kmer_spectrum[m1_range, 'frequency'],
                                                        k, kmerCoverageEst, lengthEst, rEst)

lines(predict(two_tissue_model, response = T) ~ male1_kmer_spectrum[m1_range, 'coverage'], lwd = 2, col = 'purple')
```

### BH3-2

```{R}
# MALE 2
BH3_male_k21_file = './data/genome_profiling/BH3-2/kmer_k21_full.hist'
BH3_male_k17_file = './data/genome_profiling/BH3-2/BH3-2_k17_reduced.hist'

male2_kmer_spectrum <- read.table(BH3_male_k17_file, col.names = c('coverage', 'frequency'))

m2_range <- 5:100

source('scripts/genome_models.R') # load genome models (see the script for their details)

plot(male2_kmer_spectrum[m2_range, 'frequency'] ~ male2_kmer_spectrum[m2_range, 'coverage'])

k = 17
eyeballed_coverage_m2 = 15

### original - works
# source('scripts/genome_models.R')
genomescope <- nls_4peak(male2_kmer_spectrum[m2_range, 'coverage'],
                         male2_kmer_spectrum[m2_range, 'frequency'],
                         k,
                         eyeballed_coverage_m2,
                         300e6)

biasEst <- coef(genomescope_2peak)['bias']
kmerCoverageEst<- coef(genomescope_2peak)['kmercov']
# lengths and genome size better derived from female hist
lengthEst <- mean(female_parameters[, 'length'])
rEst <- mean(female_parameters[, 'r'])

conditional_model <- nls_2peak_conditional(male2_kmer_spectrum[m2_range, 'coverage'],
                                           male2_kmer_spectrum[m2_range, 'frequency'],
                                           k,
                                           kmerCoverageEst,
                                           lengthEst,
                                           rEst,
                                           biasEst)

biasEst <- coef(genomescope_2peak)['biasEst']

m2_range <- 5:80
x <- male2_kmer_spectrum[m2_range, 'coverage']
y <- male2_kmer_spectrum[m2_range, 'frequency']
estKmercov <- 14.7666
length <- 289738623
r <- 0.003627879
biasEst <- 1
fraction_diploidEst <- coef(two_tissue_model)['fraction_diploid']
max_iterations <- 40

m2_4peak_constant_shift <- nls_4peak_two_tissue_constant_shift(x, y, 17, estKmercov, length, r, 0.7463)
lines(predict(m2_4peak_constant_shift, response = T) ~ x, lwd = 2, col = 'yellow')
nlstools::confint2(m2_4peak_constant_shift, level = 0.95, method = "asymptotic")

model2 <- nlsLM(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov + cov_shift) +
                    ((1-r)^k)        * dnbinom(x, size = (2 * kmercov) / bias, mu = (2 * kmercov) + cov_shift)) * length * fraction_diploid +
                                       dnbinom(x, size = kmercov / bias, mu = kmercov + cov_shift) * length * (1 - fraction_diploid),
                  start = list(kmercov=estKmercov, cov_shift=0, bias = biasEst, fraction_diploid=0.8),
                  lower = c("kmercov" = 0, "cov_shift" = 0, "bias" = 0, "fraction_diploid" = 0),
                  upper = c("kmercov" = Inf, "cov_shift" = Inf, "bias" = Inf, "fraction_diploid" = 1),
                  control = list(minFactor=1e-12, maxiter=max_iterations))

lines(predict(model2, response = T) ~ x)

### modeling two tissue k-mer spectra
source('scripts/genome_models.R') # load genome models (see the script for their details)

constant_shift_model <- nls_2peak_two_tissue_constant_shift(male2_kmer_spectrum[m2_range, 'coverage'],
                                                            male2_kmer_spectrum[m2_range, 'frequency'],
                                                           k,
                                                           kmerCoverageEst,
                                                           lengthEst,
                                                           rEst,
                                                           biasEst, fraction_diploidEst)
lines(predict(constant_shift_model, response = T) ~ male2_kmer_spectrum[m2_range, 'coverage'], col = 'red')
```

##########
Now, we can fix the haplid genome length (should be same in male and female) and to further reduce the number of parameters we can also assume the same heterozygosity (which won't ) really change the effect of the peak shift, which lift the burden of fit. The two new parameters will fit instead are `fraction_diploid`, which is the proportion of the genome that is autosomal and `cov_shift`, which is the coverage coresponding to the sperm tissue. The `kmercov` parameter will then have a different interpretation - the coverage of the soma tissue. We can also calculate the confidence interval for the shift.

```{R}
lengthEst <- mean(female_parameters$length)
# rEst <- mean(female_parameters$r)

nls_2peak_two_tissue_constant_shift_basic(male1_kmer_spectrum[m1_range, 'coverage'],
                                          male1_kmer_spectrum[m1_range, 'frequency'],
                                          21,
                                          50,
                                          lengthEst)

rEst <- mean(female_parameters$r)
### modeling two tissue k-mer spectra
constant_shift_model <- nls_2peak_two_tissue_constant_shift(male1_kmer_spectrum[m1_range, 'coverage'],
                                       male1_kmer_spectrum[m1_range, 'frequency'],
                                       21,
                                       50,
                                       lengthEst,
                                       rEst,
                                       biasEst)

lines(predict(constant_shift_model, response = T) ~ x, col = 'red')

nlstools::confint2(constant_shift_model, level = 0.95, method = "asymptotic")

fracEst <- coef(conditional_X0_model)['fraction_diploid']


### simplified - DOES NOT
genome_model <- nls_2peak(x, y, estKmercov, max_iterations)
lines(predict(genome_model, response = T) ~ kmer_spectrum$coverage[range])
```

## Unused exploration

let's do an R script that will mimic what we did with the mapped reads - kernel smoothing.

```
wget https://raw.githubusercontent.com/schatzlab/genomescope/master/genomescope.R -O scripts/genomescope.R
Rscript scripts/genomescope.R ./data/genome_profiling/BH3-2/kmer_k21_full.hist
```

### Kernel smoothing

I initially tried to get it work using kernel smoothing exactly in the same way as I did with mapping coverage.

```R
input_file3 = './data/genome_profiling/BH3-2/kmer_k21_full.hist'
source('scripts/get_peaks.R')

kmer_spectrum <- read.table(input_file3, col.names = c('coverage', 'frequency'))

range <- 15:150
adjust = 1

second_deriv <- diff(sign(diff(kmer_spectrum$frequency[range])))

peak_covs <- kmer_spectrum$coverage[which(second_deriv == -2) + 1]
peak_heights <- kmer_spectrum$frequency[which(second_deriv == -2) + 1]
# head(peak_covs[order(peak_heights, decreasing=T)])

ks <- density(kmer_spectrum[range, 'coverage'], bw = "SJ", adjust = adjust, weights = kmer_spectrum[range, 'frequency'] / sum(kmer_spectrum[range, 'frequency']))

plot(ks)
get_peaks(ks)
```

This approach would work for samples Ocin2 and Afus1, but it fails on BH3-2 simply because the two peaks are overlapping to much. After this exploration we focused in non-linear regression models as they allow estimating a mixture of distrubtions. In our case negative binomials.


### Basics of non-linaer least square

Just to figure how it works I generated some random data mimicking k-mer spectra and started exploring how can I fit there the distrubtions using `nls` function

```R
generated_data <- hist(rnorm(1000, 50, 10))
data <- data.frame(coverage = generated_data$mids, count = generated_data$counts)

# THIS does not work because there is no scaling term - dnorm sums only to 1 - it's a probability distribution
model1 <- nls(count ~ dnorm(coverage, mean = m, sd = s),
              start = list(m=28, s=15), data = data)

# THIS WORKS - with a normal data
fit_norm = nls(count ~ N * dnorm(coverage, m, s), data=data, start=c(m=55, s=12, N=sum(data$count)) )   

# let's generate something that is more coverage histogram like - discrete empirical distribution
generated_data <- table(round(rnorm(1000, 50, 10)))
kmer_hist <- data.frame(coverage = as.numeric(names(generated_data)), count = as.numeric(generated_data))

# yep, the normal model still works
fit_norm = nls(count ~ N * dnorm(coverage, m, s), data=kmer_hist, start=c(m=55, s=12, N=sum(kmer_hist$count)) )   

# let's try negative binomial now
fit_binom <- nls(count ~ N * dnbinom(coverage, size = size_est, mu = kmercov_1n),
									start = list(kmercov_1n=48, size_est = 50, N=sum(data$count)), data=kmer_hist)

plot(count ~ coverage, data = kmer_hist)
lines(predict(fit_binom, response = T) ~ kmer_hist$coverage)

# let's generate something bimodal, first with nearly non-overlaping peaks
generated_data <- table(round(c(rnorm(1000, 50, 10), rnorm(1000, 110, 20))))
kmer_hist <- data.frame(coverage = as.numeric(names(generated_data)), count = as.numeric(generated_data))

fit_two_binom <- nls(count ~ N * (dnbinom(coverage, size = kmercov_1n / bias, mu = kmercov_1n) + dnbinom(coverage, size = kmercov_2n / bias, mu = kmercov_2n)),
									start = list(kmercov_1n=48, kmercov_2n=105, bias = 1, N=sum(kmer_hist$count)), data=kmer_hist)


plot(count ~ coverage, data = kmer_hist)
lines(predict(fit_two_binom, response = T) ~ kmer_hist$coverage)
```

I am ready to create genome models for real data.
