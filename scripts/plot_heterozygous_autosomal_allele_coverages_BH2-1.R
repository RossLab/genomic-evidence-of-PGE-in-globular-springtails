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

pal <- c('black', rgb(0.02, 0.45, 0.65, 0.55), rgb(0.8, 0.05, 0.1, 0.55))
cex_legend <- 1.4
xlim <- c(0, 50)
ylim <- c(0, 0.1)
minor <- c(10.03, 10.03)
minor_mean <- c(11.3, 11.3)
major <- c(19.47, 19.47)
major_mean <- c(18.5, 18.5)

png('figures/het_autosomal_allele_supports/autosomal_variant_coverages_BH3-2.png')
	hist(informative_A_snps$cov_minor, col = pal[2], main = paste('Coverages supporting alleles in', ind), xlab = 'Coverage', breaks = 60, freq = F, border = F, xlim = xlim, ylim = ylim, cex.axis = cex_legend, cex.lab = cex_legend)
	hist(informative_A_snps$total_cov - informative_A_snps$cov_minor, col = pal[3], add = T, breaks = 60, freq = F, border = F)

	lines(minor, c(0, 1e6), lwd = 3, lty = 2)
	# lines(minor_mean, c(0, 1e6), lwd = 3, lty = 4)

	lines(major, c(0, 1e6), lwd = 3, lty = 2)
	# lines(major_mean, c(0, 1e6), lwd = 3, lty = 4)

	legend('topright', pch = 20, col = pal[c(NA, 2,3)], c('', 'minor (paternal) alleles', 'major (maternal) alleles'), bty = 'n', cex = cex_legend)
dev.off()

png('figures/het_autosomal_allele_supports/autosomal_and_X_variant_coverages_BH3-2.png')
	hist(informative_X_snps$total_cov, col = pal[1], main = paste('Coverages supporting alleles in', ind), xlab = 'Coverage', breaks = 100, freq = F, xlim = xlim, ylim = ylim, cex.axis = cex_legend, cex.lab = cex_legend)
	hist(informative_A_snps$cov_minor, col = pal[2], add = T, breaks = 60, freq = F, border = F)
	hist(informative_A_snps$total_cov - informative_A_snps$cov_minor, col = pal[3], add = T, breaks = 60, freq = F, border = F)

	lines(minor, c(0, 1e6), lwd = 3, lty = 2)
	lines(major, c(0, 1e6), lwd = 3, lty = 2)

	legend('topright', pch = 20, col = pal, c('X chromosome alleles', 'minor (paternal) alleles', 'major (maternal) alleles'), bty = 'n', cex = cex_legend)
dev.off()
