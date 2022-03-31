#!/usr/bin/env Rscript
ind <- 'BH3-2'
snp_tab_A <- read.table('data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_A.tsv')
snp_tab_X <- read.table('data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_X.tsv')
colnames(snp_tab_A) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')
colnames(snp_tab_X) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')

snp_tab_A$cov_minor <- apply(snp_tab_A[, c('ref_cov', 'alt_cov')], 1, min)
snp_tab_X$cov_minor <- apply(snp_tab_X[, c('ref_cov', 'alt_cov')], 1, min)

heterozygous <- sapply(snp_tab_A$genotype, function(x){gty = unlist(strsplit(x, '/')); return(gty[1] != gty[2])}  )
snp_tab_A[heterozygous, 'genotype'] <- '0/1'

ref_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/0', ]
informative_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/1' & snp_tab_A$total_cov < 100, ]
informative_X_snps <- snp_tab_X[snp_tab_X$genotype == '1/1' & snp_tab_X$total_cov < 100, ]


print(paste('A 0/0: ', nrow(snp_tab_A[snp_tab_A$genotype == '1/1' & snp_tab_A$total_cov < 100, ])))
# [1] 915879
print(paste('A 0/1: ', nrow(informative_A_snps)))
# [1] 1959258
print(paste('X 1/1: ', nrow(informative_X_snps)))
# [1] 400001
print(paste('X 0/1: ', nrow(snp_tab_X[snp_tab_X$genotype == '0/1' & snp_tab_X$total_cov < 100, ])))
# [1] 144204

figure_name = 'figures/het_autosomal_allele_supports/autosomal_and_X_variant_coverages_BH3-2.png'
source('scripts/load_palette.R')
cex_legend <- 1 #0.95
xlim <- c(0, 50)
ylim <- c(0, 0.1)

#### calculate the coverage expecation!
BH32_sperm <- read.table('tables/two_tissue_model_empirical.tsv', header = T, row.names = 1)['BH3-2', 'two_tissue_sperm']
cov_A <- mean(ref_A_snps$ref_cov)
# using the equation from the paper c_x = c_A / (2 - f_h)
cov_maternal <- cov_A / (2 - BH32_sperm)
cov_paternal <- cov_A - cov_maternal

mean_minor <- c(cov_paternal, cov_paternal) #
mean_major <- c(cov_maternal, cov_maternal) #

print(paste("Paternal expecation:", round(cov_paternal, 2)))
print(paste("Maternal expecation:", round(cov_maternal, 2)))

coverage_data <- list(X = informative_X_snps$total_cov, A_minor = informative_A_snps$cov_minor, A_major = informative_A_snps$total_cov - informative_A_snps$cov_minor)
main <- ''# paste('Coverages supporting autosomal heterozygous alleles in', ind)

source('scripts/fixed_bin_historgram.R')

plot_means = T
plot_medians = F

png(figure_name, units="in", width=5, height=5, res=300)

	# 'c(bottom, left, top, right)'
	par(mar = c(4, 4, 1, 1) + 0.1)
	fixed_bin_histogram(coverage_data, pal, main = main, xlab = 'Coverage support', bins = 50, freq = F, xlim = xlim, ylim = ylim, default_legend = F)

	lines(mean_minor, c(0, 1e6), lwd = 2, lty = 2)
	lines(mean_major, c(0, 1e6), lwd = 2, lty = 2)

	# legend('topright', pch = c(20, 20, 20, NA), lty = c(NA, NA, NA, 2), col = c(pal, 'black'), c('X chromosome alleles', 'A minor (paternal) alleles', 'A major (maternal) alleles', 'expectation'), bty = 'n', cex = cex_legend, lwd = 2)
	if ( plot_means & plot_medians){
		legend('topright', lty = c(NA, 2, 4), col = c('black'), c('', 'mean', 'median'), bty = 'n', cex = cex_legend, lwd = 2, title = 'expectations')
	} else {
		if ( plot_medians ){
			legend('topright', lty = c(NA, 4), col = c('black'), c('', 'expectation'), bty = 'n', cex = cex_legend, lwd = 2)
		} else {
			legend('topright', lty = c(NA, 2), col = c('black'), c('', 'expectation'), bty = 'n', cex = cex_legend, lwd = 2)
		}
	}

dev.off()

source('scripts/get_peaks.R')

plot_minor_allele_freq_hist <- function(minor_allele_freq, col, main, adjust = 1.5, min_cov = 0){
	ks <- density(minor_allele_freq, adjust = adjust)
	peaks <- get_peaks(ks)
	print(peaks)
	hist(minor_allele_freq, col = col, freq = F, breaks = 20, xlab = 'minor allele coverage ratio', main = main)
	lines(ks, lwd = 3)
	arrows(peaks[nrow(peaks), 'cov'], peaks[nrow(peaks), 'height'] + 0.5, peaks[nrow(peaks), 'cov'], peaks[nrow(peaks), 'height'] + 0.2, length = 0.1, lwd = 2)
	text(peaks[nrow(peaks), 'cov'], peaks[nrow(peaks), 'height'] + 0.7, round(peaks[nrow(peaks), 'cov'], 4))
}

figure_name = 'figures/het_autosomal_allele_supports/autosomal_coverage_ratio_BH3-2.png'
png(figure_name, units="in", width=5, height=5, res=300)
	plot_minor_allele_freq_hist(informative_A_snps$cov_minor / informative_A_snps$total_cov, 'grey', main = 'BH3-2', 2)
dev.off()
