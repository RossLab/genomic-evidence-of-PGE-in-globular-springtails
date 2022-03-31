source('scripts/two_tissue_model_functions.R')

######################
## A. fusca females ##
######################

females <- c('WW5-3', 'WW5-5', 'WW2-6')

WW5_3_k21_file = 'data/genome_profiling/WW5-3/kmer_k21_full.hist'
WW5_5_k21_file = 'data/genome_profiling/WW5-5/kmer_k21_full.hist'
WW2_6_k21_file = 'data/genome_profiling/WW2-6/kmer_k21_full.hist'

### female
female1_kmer_spectrum <- read.table(WW5_3_k21_file, col.names = c('coverage', 'frequency'))
female2_kmer_spectrum <- read.table(WW5_5_k21_file, col.names = c('coverage', 'frequency'))
female3_kmer_spectrum <- read.table(WW2_6_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:120
# plot(frequency ~ coverage, data = female2_kmer_spectrum[cov_range, ])

#
female1_model <- nls_4peak(female1_kmer_spectrum[cov_range, 'coverage'], female1_kmer_spectrum[cov_range, 'frequency'], 21, 15, 300e6, 40)
female2_model <- nls_4peak(female2_kmer_spectrum[cov_range, 'coverage'], female2_kmer_spectrum[cov_range, 'frequency'], 21, 20, 300e6, 40)
female3_model <- nls_4peak(female3_kmer_spectrum[cov_range, 'coverage'], female3_kmer_spectrum[cov_range, 'frequency'], 21, 15, 300e6, 40)

female_parameters <- as.data.frame(rbind(coef(female1_model), coef(female2_model), coef(female3_model)))
female_parameters$ind <- females
rownames(female_parameters) <- females

write.table(female_parameters, 'tables/kmer_spectra_modeling_female_parameters.tsv', quote = F, col.names = T, row.names = F)
# these are practically genomescope estimates, we can also fit there some of the models that fit indipendently 1n and 2n coverages

# FITTING FEMALE PGE MODELS
# # here, using regular GenomeScope estimates of length and heterozygosity as fixed parameters to reduce complexity of the models
female1_PGE_model <- nlsLM_2peak_proportional_peaks(female1_kmer_spectrum[cov_range, 'coverage'], female1_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)

female2_PGE_model <- nlsLM_2peak_proportional_peaks(female2_kmer_spectrum[cov_range, 'coverage'], female2_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)

cov_range <- 7:120 # this one has a bit wider error profile
female3_PGE_model <- nlsLM_2peak_proportional_peaks(female3_kmer_spectrum[cov_range, 'coverage'], female3_kmer_spectrum[cov_range, 'frequency'], 15, 280e6, 0.2)



############################
## Orchesella cincta male ##
############################

Ocin2_k21_file = './data/genome_profiling/Ocin2/kmer_k21_full.hist'
Ocin2_kmer_spectrum <- read.table(Ocin2_k21_file, col.names = c('coverage', 'frequency'))
cov_range <- 5:150

Ocin2_Genomescope <- nls_4peak(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 21, 50, 280e6, 40)

Ocin2_PGE_model <- nlsLM_2peak_proportional_peaks(Ocin2_kmer_spectrum[cov_range, 'coverage'],  Ocin2_kmer_spectrum[cov_range, 'frequency'], 50, 300e6, 0.3)

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
# source('scripts/genome_models.R')
Afus1_genomescope <- nls_4peak(Afus1_kmer_spectrum[Afus1_range, 'coverage'],
                         Afus1_kmer_spectrum[Afus1_range, 'frequency'],
                         k, eyeballed_Afus1_coverage, 300e6)

Afus1_PGE_model <- nlsLM_2peak_proportional_peaks(Afus1_kmer_spectrum[Afus1_range, 'coverage'],  Afus1_kmer_spectrum[Afus1_range, 'frequency'], 55, 280e6, 0.5)

# MALE Bh3-2
# BH3_male_k21_file = './data/genome_profiling/BH3-2/kmer_k21_full.hist'
BH3_male_k17_file = './data/genome_profiling/BH3-2/BH3-2_k17_reduced.hist'

BH32_kmer_spectrum <- read.table(BH3_male_k17_file, col.names = c('coverage', 'frequency'))
BH32_range <- 5:100

k = 17
eyeballed_coverage_BH32 = 15

### original - works
# source('scripts/genome_models.R')
BH32_genomescope <- nls_4peak(BH32_kmer_spectrum[BH32_range, 'coverage'],
                              BH32_kmer_spectrum[BH32_range, 'frequency'],
                              k,
                              eyeballed_coverage_BH32,
                              300e6)

BH32_PGE_model <- nlsLM_2peak_proportional_peaks(BH32_kmer_spectrum[BH32_range, 'coverage'],  BH32_kmer_spectrum[BH32_range, 'frequency'], eyeballed_coverage_BH32, 280e6, 0.5)
###################
## MODEL SUMMARY ##
###################

inds <- c('Afus1', 'BH3-2', 'Ocin2', females)
GenomeScope_models <- list(Afus1_genomescope, BH32_genomescope, Ocin2_Genomescope, female1_model, female2_model, female3_model)
PGE_models <- list(Afus1_PGE_model, BH32_PGE_model, Ocin2_PGE_model, female1_PGE_model, female2_PGE_model, female3_PGE_model)

PGE_2n_coverages <- sapply(PGE_models, function(m){ coef(m)['kmercov2'] })
PGE_1n_coverages <- PGE_2n_coverages * sapply(PGE_models, function(m){ coef(m)['proportion'] })



two_tissure_model_frame <- data.frame(ind = inds, 'cov_1n' = round(PGE_1n_coverages, 2), 'cov_2n' = round(PGE_2n_coverages, 2))
# save the estimates?

two_tissure_model_frame$two_tissue_sperm <- round(1 - ((PGE_2n_coverages - PGE_1n_coverages) / PGE_1n_coverages), 3)
two_tissure_model_frame[, 'overdispersal'] <- round(sapply(PGE_models, function(m){ coef(m)['bias'] }), 3)

# confidence intervals?
# nlstools::confint2(Afus1_PGE_model, level = 0.95, method = "asymptotic")
# nlstools::confint2(BH32_PGE_model, level = 0.95, method = "asymptotic")
# nlstools::confint2(Ocin2_PGE_model, level = 0.95, method = "asymptotic")

CIs <- t(sapply(PGE_models, nlstools::confint2, level = 0.95, method = "asymptotic"))
two_tissure_model_frame[, 'cov_1n/cov_2n'] <- round(sapply(PGE_models, function(m){ coef(m)['proportion'] }), 3)
two_tissure_model_frame[, 'cov_1n/cov_2n_CI_l'] <- round(CIs[, 1], 3)
two_tissure_model_frame[, 'cov_1n/cov_2n_CI_u'] <- round(CIs[, 6], 3)

write.table(two_tissure_model_frame, 'tables/two_tissue_model_empirical.tsv', col.names = T, quote = F, row.names = F, sep = '\t')

###################
## VISUALISATION ##
###################

pdf('figures/two_tissue_model_BH3-2.pdf')
  plot_PGE_model(BH32_PGE_model, xlim = c(5, 70), F, 1.3, bty = 'n', legend = F)
dev.off()

# plot(Afus1_kmer_spectrum[Afus1_range, 'frequency'] ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'])
# lines(predict(Afus1_genomescope, response = T) ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'])
# lines(predict(Afus1_PGE_model, response = T) ~ Afus1_kmer_spectrum[Afus1_range, 'coverage'], col = 'red', lty = 2, pwd = 2)

pdf('figures/two_tissue_model_BH3-2_SM.pdf')
  plot_PGE_model(BH32_PGE_model, xlim = c(0, 70), T, 1.3, bty = 'n', legend = F)
dev.off()

pdf('figures/two_tissue_model_Afus1_SM.pdf')
  plot_PGE_model(Afus1_PGE_model, xlim = c(0, 220), T, 1.3, bty = 'n', legend = F)
dev.off()

pdf('figures/two_tissue_model_Ocin2_SM.pdf')
  plot_PGE_model(Ocin2_PGE_model, xlim = c(0, 150), T, 1.3, bty = 'n', legend = T)
dev.off()

# females
# pdf('figures/two_tissue_model_WW5-3.pdf')
#   plot_PGE_model(female1_PGE_model)
# dev.off()
#
# pdf('figures/two_tissue_model_WW5-5.pdf')
#   plot_PGE_model(female2_PGE_model)
# dev.off()
#
# pdf('figures/two_tissue_model_WW2-6.pdf')
#   plot_PGE_model(female3_PGE_model)
# dev.off()
