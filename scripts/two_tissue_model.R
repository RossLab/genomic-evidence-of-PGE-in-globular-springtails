#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
#     help="Print extra output [default]")
parser$add_argument("-i", "--input", dest="kmer_hist_file",
                    help="Input histogram file")

args <- parser$parse_args()

source('scripts/two_tissue_model_functions.R') # functions

# kmer_spectrum <- read.table('data/generated/kmcdb_k21.hist', col.names = c('coverage', 'frequency'))
kmer_spectrum <- read.table(args$kmer_hist_file, col.names = c('coverage', 'frequency'))

cov_range <- 5:120

GenomeScope_model <- nls_4peak(kmer_spectrum[cov_range, 'coverage'], kmer_spectrum[cov_range, 'frequency'], 21, 28, 20e6, 40)
PGE_model <- nlsLM_2peak_unconditional_peaks(kmer_spectrum[cov_range, 'coverage'], kmer_spectrum[cov_range, 'frequency'], 25, 10e6, 0.2)

plot(kmer_spectrum[cov_range, 1], kmer_spectrum[cov_range, 2])
