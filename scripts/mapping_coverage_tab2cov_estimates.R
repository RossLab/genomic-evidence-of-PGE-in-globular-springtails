#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
#     help="Print extra output [default]")
parser$add_argument("-i", "--input", dest="coverage_tab_file",
                    help="Input coverage table file")
parser$add_argument("-o", "--output_pattern", dest="out",
                    help="The output pattern for the summary and plot")

args <- parser$parse_args()

coverage <- read.table(args$coverage_tab_file)

scf_ks <- density(coverage[, 3], bw = "SJ", adjust = 1)
second_deriv <- diff(sign(diff(scf_ks$y)))

peak_covs <- scf_ks$x[which(second_deriv == -2) + 1]
peak_heights <- scf_ks$y[which(second_deriv == -2) + 1]


# male
est_1 <- peak_covs[which.max(peak_heights)]
est_2 <- peak_covs[peak_heights!=max(peak_heights)][which.max( peak_heights[peak_heights!=max(peak_heights)] )]
if (est_1 > est_2){
		monosomic_coverage <- est_2
		disomic_coverage <- est_1
} else {
		monosomic_coverage <- est_1
		disomic_coverage <- est_2
}

output_summary_file <- paste0(args$out, '_summary_tsv')
write(paste(c("mapping_monosomic", "mapping_disomic"), collapse = '\t'), output_summary_file)
write(paste(c(monosomic_coverage, disomic_coverage), collapse = '\t'), output_summary_file, append = T)

output_plot <- paste0(args$out, '_plot.png')
png(output_plot)
  plot(scf_ks)
dev.off()
