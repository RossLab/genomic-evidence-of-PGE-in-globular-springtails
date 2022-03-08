#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
#     help="Print extra output [default]")
parser$add_argument("-i", "--input", dest="kmer_hist_file",
                    help="Input histogram file")
parser$add_argument("-o", "--output_pattern", dest="out",
                    help="The output pattern for the summary and plot")

args <- parser$parse_args()

source('scripts/two_tissue_model_functions.R') # functions

# kmer_spectrum <- read.table('data/generated/kmcdb_k21.hist', col.names = c('coverage', 'frequency'))
kmer_spectrum <- read.table(args$kmer_hist_file, col.names = c('coverage', 'frequency'))

write("loaded k-mer spectra", stderr())

# removing errors from the model (removing all monotonicly decreasing values)
min_coverage <- min(which(diff(kmer_spectrum$frequency) > 0)) - 1
# and capping 99.9% of kmers in the dataset excluding erros
max_coverage <- min(which(cumsum(kmer_spectrum$frequency[min_coverage:nrow(kmer_spectrum)]) / sum(kmer_spectrum$frequency[min_coverage:nrow(kmer_spectrum)]) > 0.999))
cov_range <- min_coverage:max_coverage


# to get coverage mode I just use second derivation
second_deriv <- diff(sign(diff(kmer_spectrum$frequency[cov_range])))
peak_covs <- kmer_spectrum$coverage[cov_range][which(second_deriv == -2) + 1]
peak_heights <- kmer_spectrum$frequency[cov_range][which(second_deriv == -2) + 1]

coverage_mode <- peak_covs[which.max(peak_heights)]
coverage_mean <- weighted.mean(kmer_spectrum$coverage[cov_range], kmer_spectrum$frequency[cov_range])

if (coverage_mode < coverage_mean){
  cov_prior <- coverage_mode
} else {
  cov_prior <- (coverage_mode / 2)
}

length_prior <- sum(kmer_spectrum$coverage[cov_range] * kmer_spectrum$frequency[cov_range]) / cov_prior

print("Priors: ")
print(c(min_coverage, max_coverage))
print(cov_prior)
print(length_prior)

GenomeScope_model <- nls_4peak(kmer_spectrum[cov_range, 'coverage'], kmer_spectrum[cov_range, 'frequency'], 21, cov_prior, length_prior)
PGE_model <- nlsLM_2peak_unconditional_peaks(kmer_spectrum[cov_range, 'coverage'], kmer_spectrum[cov_range, 'frequency'], cov_prior, length_prior, 0.5)

###
# summary of results
###

genomescope_estimates <- coef(GenomeScope_model)
PGE_estimates <- coef(PGE_model)

write("GenomeScope estimates: ", stderr())
write(paste(names(genomescope_estimates), '=', round(genomescope_estimates, 4)), stderr())
write("Two tissue model estimates: ", stderr())
write(paste(names(PGE_estimates), '=', round(PGE_estimates, 4)), stderr())

write(paste(c(genomescope_estimates, PGE_estimates), sep = '\t'), stdout())

output_summary_file <- paste0(args$out, '_summary_tsv')
write(paste(c(names(genomescope_estimates), names(PGE_estimates)), collapse = '\t'), output_summary_file)
write(paste(c(genomescope_estimates, PGE_estimates), collapse = '\t'), output_summary_file, append = T)

###
# plot
###

output_plot <- paste0(args$out, '_plot.png')
png(output_plot)
  plot(kmer_spectrum[cov_range, 1], kmer_spectrum[cov_range, 2], xlab = 'Coverage', ylab = 'Frequency')

  lines(kmer_spectrum[cov_range, 1], predict(GenomeScope_model), lty = 1, lwd = 2, col = 'blue')
  lines(kmer_spectrum[cov_range, 1], predict(PGE_model), lty = 2, lwd = 2, col = 'purple')

  legend('topright', c('simulated k-mer spectra', 'GenomeScope model', 'Two tissue model'), pch = c(1, NA, NA), lty = c(NA, 1, 2), lwd = 2, col = c('black', 'blue', 'purple'), bty = 'n')
dev.off()
