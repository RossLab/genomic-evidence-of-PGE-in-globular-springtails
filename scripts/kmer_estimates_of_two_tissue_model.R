require('minpack.lm')

###############
## FUNCTIONS ##
###############

### Original genomescope model
nls_4peak <- function(x, y, k, estKmercov, estLength, max_iterations = 40){
    model4 = NULL

    cat("Original genomescope model 'nls_4peak'\n")

    try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                          (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                          (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                          (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                      start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)

    if(class(model4) == "try-error"){

        try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                              (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                              (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                              (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                          start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                          algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)
    }

    return(model4)
}

# two tissue model
nlsLM_2peak_unconditional_peaks <- function(x, y, kmerEst, lengthEst, hetEst = 0.6){
  nlsLM(y ~ ((het       * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
            ((1 - het)  * dnbinom(x, size = kmercov2  / bias, mu = kmercov2))) * length,
        start = list(kmercov = kmerEst, kmercov2 = (2 * kmerEst), bias = 0.5, length = lengthEst, het = hetEst),
        control = list(minFactor=1e-12, maxiter=40))
}

predict_1n_peak_nlsLM_2peak_unconditional <- function(model){
  model_env <- model$m$getEnv()
  model_env$het * dnbinom(model_env$x, size = model_env$kmercov   / model_env$bias, mu = model_env$kmercov) * model_env$length
}

predict_2n_peak_nlsLM_2peak_unconditional <- function(model){
  model_env <- model$m$getEnv()
  (1 - model_env$het)  * dnbinom(model_env$x, size = model_env$kmercov2  / model_env$bias, mu = model_env$kmercov2) * model_env$length
}

coverage_barplot <- function(bar_heights, bar_positions, font_size = 1, width = 0.5){

  plot(bar_heights, type="n", xlab="Coverage", ylab="Frequency",
       ylim=c(0, max(bar_heights)), xlim=range(bar_positions),
       cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  for ( i in 1:length(bar_heights)){
    rect(bar_positions[i] - width, 0, bar_positions[i] + width, bar_heights[i], col = 'deepskyblue', border = F)
  }
}

plot_PGE_model <- function(model){
  x <- model$m$getEnv()$x
  y <- model$m$getEnv()$y
  cov_1n <- coef(model)['kmercov']
  cov_2n <- coef(model)['kmercov2']
  monoploid_col <- 'darkorchid4'
  diploid_col <- 'chocolate'

  coverage_barplot(y, x)

  disomic_prediction <- predict_2n_peak_nlsLM_2peak_unconditional(model)
  monosomic_prediction <- predict_1n_peak_nlsLM_2peak_unconditional(model)

  lines(predict(model, response = T) ~ x, lwd = 3)
  lines(monosomic_prediction ~ x, lwd = 3, col = monoploid_col)
  lines(disomic_prediction ~ x, lwd = 3, col = diploid_col)

  lines(c(cov_1n, cov_1n) - 0.5, c(0, max(y)), lty = 2)
  lines(c(cov_2n, cov_2n) - 0.5, c(0, max(y)), lty = 2)

  legend('topright',
         c('kmer histogram','full model', 'monoploid', 'diploid'),
         col = c('deepskyblue','black', monoploid_col, diploid_col),
         lty = c(NA, 1, 1, 1), lwd = 3,
         pch = c(15, NA, NA, NA), bty = 'n')
}



######################
## A. fusca females ##
######################

females <- c('WW5-3', 'WW5-5', 'WW2-6')

WW5_3_k21_file = './data/genome_profiling/WW5-3/kmer_k21_full.hist'
WW5_5_k21_file = './data/genome_profiling/WW5-5/kmer_k21_full.hist'
WW2_6_k21_file = './data/genome_profiling/WW2-6/kmer_k21_full.hist'

### female
feAfus1_kmer_spectrum <- read.table(WW5_3_k21_file, col.names = c('coverage', 'frequency'))
feBH32_kmer_spectrum <- read.table(WW5_5_k21_file, col.names = c('coverage', 'frequency'))
female3_kmer_spectrum <- read.table(WW2_6_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:120
# plot(frequency ~ coverage, data = feBH32_kmer_spectrum[cov_range, ])

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

cov_range <- 7:120 # this one has a bit wider error profile
female3_PGE_model <- nlsLM_2peak_unconditional_peaks(female3_kmer_spectrum[cov_range, 'coverage'], female3_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)



############################
## Orchesella cincta male ##
############################

Ocin2_k21_file = './data/genome_profiling/Ocin2/kmer_k21_full.hist'
Ocin2_kmer_spectrum <- read.table(Ocin2_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:150

Ocin2_Genomescope <- nls_4peak(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 21, 50, 280e6, 40)

Ocin2_PGE_model <- nlsLM_2peak_unconditional_peaks(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 50, 300e6, 0.3)

####################
## A. fusca males ##
####################

# REFERENCE MALE Afus1
Afus1_k21_file = './data/genome_profiling/Afus1/kmer_k21_full.hist'
Afus1_kmer_spectrum <- read.table(Afus1_k21_file, col.names = c('coverage', 'frequency'))
Afus1_range <- 20:250

k = 21
eyeballed_Afus1_coverage = 50

### original - works
Afus1_genomescope <- nls_4peak(Afus1_kmer_spectrum[Afus1_range, 'coverage'],
                         Afus1_kmer_spectrum[Afus1_range, 'frequency'],
                         k, eyeballed_Afus1_coverage, 300e6)

Afus1_PGE_model <- nlsLM_2peak_unconditional_peaks(Afus1_kmer_spectrum[Afus1_range, 'coverage'],  Afus1_kmer_spectrum[Afus1_range, 'frequency'], 55, 280e6, 0.5)

# MALE Bh3-2
# BH3_male_k21_file = './data/genome_profiling/BH3-2/kmer_k21_full.hist'
BH3_male_k17_file = './data/genome_profiling/BH3-2/BH3-2_k17_reduced.hist'

BH32_kmer_spectrum <- read.table(BH3_male_k17_file, col.names = c('coverage', 'frequency'))
BH32_range <- 5:100

k = 17
eyeballed_coverage_BH32 = 15

### original - works
BH32_genomescope <- nls_4peak(BH32_kmer_spectrum[BH32_range, 'coverage'],
                              BH32_kmer_spectrum[BH32_range, 'frequency'],
                              k,
                              eyeballed_coverage_BH32,
                              300e6)

BH32_PGE_model <- nlsLM_2peak_unconditional_peaks(BH32_kmer_spectrum[BH32_range, 'coverage'],  BH32_kmer_spectrum[BH32_range, 'frequency'], eyeballed_coverage_BH32, 280e6, 0.5)

###################
## MODEL SUMMARY ##
###################

inds <- c('Afus1', 'BH3-2', 'Ocin2', females)
GenomeScope_models <- list(Afus1_genomescope, BH32_genomescope, Ocin2_Genomescope, female1_model, female2_model, female3_model)
PGE_models <- list(Afus1_PGE_model, BH32_PGE_model, Ocin2_PGE_model, female1_PGE_model, female2_PGE_model, female3_PGE_model)

PGE_1n_coverages <- sapply(PGE_models, function(m){ coef(m)['kmercov'] })
PGE_2n_coverages <- sapply(PGE_models, function(m){ coef(m)['kmercov2'] })

data.frame(ind = inds, 'cov_1n' = PGE_1n_coverages, 'cov_2n' = PGE_2n_coverages)
# save the estimates?

# confidence intervals?
# nlstools::confint2(BH32_PGE_model, level = 0.95, method = "asymptotic")

###################
## VISUALISATION ##
###################

# plot(Afus1_kmer_spectrum[Afus1_range, 'frequency'] ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'])
# lines(predict(Afus1_genomescope, response = T) ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'])
# lines(predict(Afus1_PGE_model, response = T) ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'], col = 'red', lty = 2, pwd = 2)

png('figures/Ocin2_PGE_model.png')
	plot_PGE_model(Ocin2_PGE_model)
  lines(predict(Ocin2_Genomescope, response = T) ~ Ocin2_Genomescope$m$getEnv()$x, lwd = 3, col = 'grey', lty = 3)
dev.off()

png('figures/Afus1_PGE_model.png')
	plot_PGE_model(Afus1_PGE_model)
  lines(predict(Afus1_genomescope, response = T) ~ Afus1_genomescope$m$getEnv()$x, lwd = 3, col = 'grey', lty = 3)
dev.off()

png('figures/BH3-2_PGE_model.png')
	plot_PGE_model(BH32_PGE_model)
dev.off()

png('figures/WW5-3_PGE_model.png')
	plot_PGE_model(female1_PGE_model)
dev.off()

png('figures/WW5-5_PGE_model.png')
	plot_PGE_model(female2_PGE_model)
dev.off()

png('figures/WW2-6_PGE_model.png')
	plot_PGE_model(female3_PGE_model)
dev.off()
