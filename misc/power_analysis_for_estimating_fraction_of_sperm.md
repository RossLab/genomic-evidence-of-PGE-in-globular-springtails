## redoing BH3-2 kmer spectra with k = 17

Lower k helps with low coverage datasets in cost of resolution. With k = 17, a smaller fraction of genome occurs on unique k-mers, but the k-mer coverage is (R - k + 1 / R) * C, hence lower k is, higher coverage we get. The improvement is only marginal, but we needs only a marginal improvement to make the fit easier.

```bash
conda activate ntcard
# gospring
qsub -o logs -e logs -cwd -N ntcard -V -pe smp64 16 -b yes 'time ntcard -t 16 -k 17 -c 10000 -o data/resequencing/genome_profiling/BH3-2/BH3-2_k17.hist @data/resequencing/genome_profiling/BH3-2/FILES'
```

## k-mer spectra to 1n and 2n estimates

GenomeScope fits a model that is... unsuitable for estimating 1n and 2n coverage indipendently.

Testing data

```
./data/genome_profiling/Ocin2/kmer_k21_full.hist
./data/genome_profiling/WW5-3/kmer_k21_full.hist
./data/genome_profiling/BH3-2/kmer_k21_full.hist
./data/genome_profiling/BH3-2/BH3-2_k17_reduced.hist
```

let's do an R script that will mimic what we did with the mapped reads - kernel smoothing.

```
wget https://raw.githubusercontent.com/schatzlab/genomescope/master/genomescope.R -O scripts/genomescope.R
Rscript scripts/genomescope.R ./data/genome_profiling/BH3-2/kmer_k21_full.hist
```

### Kernel smoothing

```R
input_file1 = './data/genome_profiling/Ocin2/kmer_k21_full.hist'
input_file2 = './data/genome_profiling/WW5-3/kmer_k21_full.hist'
input_file3 = './data/genome_profiling/BH3-2/kmer_k21_full.hist'
BH3_male_k17_file = './data/genome_profiling/BH3-2/BH3-2_k17_reduced.hist'


source('scripts/get_peaks.R')

kmer_spectrum <- read.table(input_file1, col.names = c('coverage', 'frequency'))

range <- 15:150
adjust = 1

second_deriv <- diff(sign(diff(kmer_spectrum$frequency[range])))

peak_covs <- kmer_spectrum$coverage[which(second_deriv == -2) + 1]
peak_heights <- kmer_spectrum$frequency[which(second_deriv == -2) + 1]
# head(peak_covs[order(peak_heights, decreasing=T)])

ks <- density(kmer_spectrum[range, 'coverage'], bw = "SJ", adjust = adjust, weights = kmer_spectrum[range, 'frequency'] / sum(kmer_spectrum[range, 'frequency']))
``
# plot(ks)
# get_peaks(ks)
# kmer_spectrum$frequency[range]
# kmeans(, centers=c(15,30))
```

### nls

```{R}
source('scripts/genome_models.R')

### female
female_kmer_spectrum <- read.table(input_file2, col.names = c('coverage', 'frequency'))
cov_range <- 5:120
plot(frequency ~ coverage, data = female_kmer_spectrum[cov_range, ])

female_model <- nls_4peak(female_kmer_spectrum[cov_range, 'coverage'], female_kmer_spectrum[cov_range, 'frequency'], 21, 15, 350e6, 40)
female_parameters <- coef(female_model)

### male
# BH3-2 k=21 input_file3
# BH3-2 k=17 BH3_male_k17_file
kmer_spectrum <- read.table(BH3_male_k17_file, col.names = c('coverage', 'frequency'))
range <- 5:60

# plot(kmer_spectrum[range, 'frequency'] * kmer_spectrum[range, 'coverage'])
plot(frequency ~ coverage, data = kmer_spectrum[range, ])

x <- kmer_spectrum[range, 'coverage']
y <- kmer_spectrum[range, 'frequency']
estKmercov <- 15
max_iterations <- 40
VERBOSE=1
k = 17
### original - works
genomescope <- nls_4peak(x, y, k, estKmercov, sum(y * x) / estKmercov, max_iterations)

### a bit simplified - WORKS
genomescope_2peak <- nls_2peak_kmer_explicit(x, y, k, estKmercov, sum(y * x) / estKmercov, max_iterations)
lines(predict(genomescope_2peak, response = T) ~ kmer_spectrum$coverage[range])

biasEst <- coef(genomescope_2peak)['bias']
kmerCoverageEst<- coef(genomescope_2peak)['kmercov']
# lengths and genome size better derived from female hist
lengthEst <- female_parameters['length']
rEst <- female_parameters['r']

conditional_model <- nls_2peak_conditional(x, y, k, kmerCoverageEst, lengthEst, rEst, biasEst, max_iterations)
lines(predict(conditional_model, response = T) ~ kmer_spectrum$coverage[range], col = 'green')
# THis conditional model is a tiny bit strange - it does not really acknowledge there is a haploid portion of genome.

conditional_X0_model <- nls_2peak_male_X0_conditional(x, y, k, kmerCoverageEst, lengthEst, rEst, biasEst, max_iterations)
# this model should be a bit more explicit
lines(predict(conditional_X0_model, response = T) ~ x, col = 'yellow')

### modeling two tissue k-mer spectra
constant_shift_model <- nls_2peak_two_tissue_constant_shift(x, y, k, kmerCoverageEst, lengthEst, rEst, biasEst, max_iterations)
lines(predict(constant_shift_model, response = T) ~ x, col = 'red')

nlstools::confint2(constant_shift_model, level = 0.95, method = "asymptotic")

fracEst <- coef(conditional_X0_model)['fraction_diploid']

### finally, the perfect-most model - explicitely models PGE in the k-mer spectra
PGE_model <- nls_3peak_two_tissue_PGE(x, y, k, estKmercov, lengthEst, rEst, biasEst, fracEst, max_iterations)
PGE_parameters <- coef(PGE_model)

PGE_parameters['kmercov_m']
PGE_parameters['kmercov_p']

lines(predict(PGE_model, response = T) ~ x, col = 'purple', lty = 2, lwd = 3)

### simplified - DOES NOT
genome_model <- nls_2peak(x, y, estKmercov, max_iterations)
lines(predict(genome_model, response = T) ~ kmer_spectrum$coverage[range])

# confint2(genome_model, level = 0.95, method = c("asymptotic"))
# est_1n <- coef(genome_model)[1]
# est_2n <- coef(genome_model)[2]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # 1n / 2n kmer coverage estimates to genome coverage
# k = 21
# readlength = 150
# kmer2genome_coverage <- function(est, k, readlength){
# 	return((est * readlength) / (readlength - k + 1))
# }
#
# kmer2genome_coverage(est_1n, k, readlength)
# kmer2genome_coverage(est_2n, k, readlength)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

########
# Quest to find reasonable starting values
########

## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
numofReads = sum(as.numeric(x*y))/(readlength-k+1)
estKmercov  = x[which(y==max(y))][1]
estCoverage1 = estKmercov1*readlength/(readlength-k)
estLength   = numofReads*readlength/estCoverage1
max_iterations = 4
VERBOSE = T
nls1    = nls_4peak(x, y, k, estKmercov1, estLength1, max_iterations)

## Second we half the max kmercoverage (typically the heterozygous peak)
estKmercov2  = estKmercov1 / 2 ##2.5
estCoverage2 = estKmercov2*readlength/(readlength-k)
estLength2   = numofReads*readlength/estCoverage2

nls2 = nls_4peak(x, y, k, estKmercov2, estLength2, max_iterations)


```




### Restart

```R
generated_data <- hist(rnorm(1000, 50, 10))
data <- data.frame(coverage = generated_data$mids, count = generated_data$counts)

# THIS does not work because there is no scaling term - dnorm won't ever generate value over
model1 <- nls(count ~ dnorm(coverage, mean = m, sd = s),
              start = list(m=28, s=15), data = data)

# THIS WORKS - with a normal data
fit_norm = nls(count ~ N * dnorm(coverage, m, s), data=data, start=c(m=55, s=12, N=sum(data$count)) )   

# let's generate something that is more coverage histogram like - discrete empirical distribution
generated_data <- table(round(rnorm(1000, 50, 10)))
kmer_hist <- data.frame(coverage = as.numeric(names(generated_data)), count = as.numeric(generated_data))

# yep, the normal model still works
fit_norm = nls(count ~ N * dnorm(coverage, m, s), data=kmer_hist, start=c(m=55, s=12, N=sum(kmer_hist$count)) )   

# this GOES THROUGH
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

confint2(fit_two_binom, level = 0.95, method = c("asymptotic"))

source('scripts/genomescope_simplified.R')

nls_2peak(kmer_hist$coverage, kmer_hist$count, 50, 50)
```
