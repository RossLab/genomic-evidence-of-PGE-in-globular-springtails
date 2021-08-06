#!/usr/bin/env Rscript
ind <- 'BH3-2'
snp_tab_A <- read.table('data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_A.tsv')
snp_tab_X <- read.table('data/SNP_calls/freebayes_Afus_filt_sorted_BH3-2_X.tsv')
colnames(snp_tab_A) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')
colnames(snp_tab_X) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')

snp_tab_A$cov_minor <- apply(snp_tab_A[, c('ref_cov', 'alt_cov')], 1, min)
snp_tab_X$cov_minor <- apply(snp_tab_X[, c('ref_cov', 'alt_cov')], 1, min)

ref_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/0', ]
informative_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/1' & snp_tab_A$total_cov < 100, ]
informative_X_snps <- snp_tab_X[snp_tab_X$genotype == '0/0' & snp_tab_X$total_cov < 100, ]

figure_name = 'figures/het_autosomal_allele_supports/autosomal_and_X_variant_coverages_BH3-2.png'
pal <- c(rgb(0, 0, 0, 0.7), rgb(0.8, 0.05, 0.1, 0.55), rgb(0.02, 0.45, 0.65, 0.55))
cex_legend <- 0.75
xlim <- c(0, 50)
ylim <- c(0, 0.1)
minor <- c(10.03, 10.03)
major <- c(19.47, 19.47)

coverage_data <- list(X = informative_X_snps$total_cov, A_minor = informative_A_snps$cov_minor, A_major = informative_A_snps$total_cov - informative_A_snps$cov_minor)
main <- ''# paste('Coverages supporting autosomal heterozygous alleles in', ind)

source('scripts/fixed_bin_historgram.R')

png(figure_name, units="in", width=5, height=5, res=300)

	# 'c(bottom, left, top, right)'
	par(mar = c(4, 4, 1, 1) + 0.1)
	fixed_bin_histogram(coverage_data, pal, main = main, xlab = 'Coverage support', bins = 50, freq = F, xlim = xlim, ylim = ylim, default_legend = F)

	lines(minor, c(0, 1e6), lwd = 3, lty = 2)
	lines(major, c(0, 1e6), lwd = 3, lty = 2)

	legend('topright', pch = 20, col = pal, c('X chromosome alleles', 'A minor (paternal) alleles', 'A major (maternal) alleles'), bty = 'n', cex = cex_legend)
dev.off()
