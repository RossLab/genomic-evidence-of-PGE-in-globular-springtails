# coverages <- read.table('data/mapped_reads/per_scf_cov_medians_Ocin2.tsv', col.names = c('scf', 'cov'))
# hist(coverages$cov[coverages$cov < 150], breaks = 150)
# 30 - 130 is reasonable cov range

ind <- 'Ocin2'
snp_tab_A <- read.table('data/SNP_calls/freebayes_Ocin2_filt_sorted_A_Ocin2.tsv')
snp_tab_X <- read.table('data/SNP_calls/freebayes_Ocin2_filt_sorted_X_Ocin2.tsv')
colnames(snp_tab_A) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')
colnames(snp_tab_X) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')


snp_tab_A$cov_minor <- apply(snp_tab_A[, c('ref_cov', 'alt_cov')], 1, min)
snp_tab_X$cov_minor <- apply(snp_tab_X[, c('ref_cov', 'alt_cov')], 1, min)

informative_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/1' & snp_tab_A$total_cov < 120 & snp_tab_A$total_cov > 50, ]
informative_X_snps <- snp_tab_X[snp_tab_X$genotype == '1/1' & snp_tab_X$total_cov < 120 & snp_tab_X$total_cov > 5, ]
coverage_data <- list(X = informative_X_snps$total_cov, A_minor = informative_A_snps$cov_minor, A_major = informative_A_snps$total_cov - informative_A_snps$cov_minor)

print(paste('A 0/0: ', nrow(snp_tab_A[snp_tab_A$genotype == '1/1' & snp_tab_A$total_cov < 120 & snp_tab_A$total_cov > 50, ])))
# [1] 915879
print(paste('A 0/1: ', nrow(informative_A_snps)))
# [1] 1959258
print(paste('X 1/1: ', nrow(informative_X_snps)))
# [1] 400001
print(paste('X 0/1: ', nrow(snp_tab_X[snp_tab_X$genotype == '0/1' & snp_tab_X$total_cov < 120 & snp_tab_X$total_cov > 50, ])))
# [1] 144204


# random_subset <- sample(1:nrow(informative_A_snps), 50000)
# library(hexbin)
# library(RColorBrewer)
# # plot(informative_A_snps$total_cov[random_subset] ~ informative_A_snps$cov_minor[random_subset])
#
# # Make the plot
# bin <- hexbin(informative_A_snps$total_cov[random_subset], informative_A_snps$cov_minor[random_subset], xbins=40)
# my_colors <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
# plot(bin, main="" , colramp=my_colors , legend=F )

figure_name = 'figures/het_autosomal_allele_supports/autosomal_and_X_variant_coverages_Ocin2.png'
source('scripts/load_palette.R')
cex_legend <- 0.95
xlim <- c(0, 120)
# expectation <- c(19.47, 19.47)
main <- ''# paste('Coverages supporting autosomal heterozygous alleles in', ind)

source('scripts/fixed_bin_historgram.R')

png(figure_name, units="in", width=5, height=5, res=300)

# 'c(bottom, left, top, right)'
par(mar = c(4, 4, 1, 1) + 0.1)
fixed_bin_histogram(coverage_data, pal, main = main, xlab = 'Coverage support', xlim = xlim, bins = 50, freq = F, default_legend = F)

# lines(expectation, c(0, 1e6), lwd = 3, lty = 2)
# legend('topright', pch = 20, col = pal, c('X chromosome alleles', 'A minor alleles', 'A major alleles'), bty = 'n', cex = cex_legend)

dev.off()

# hist(informative_A_snps$cov_minor / informative_A_snps$total_cov)

source('scripts/get_peaks.R')

plot_minor_allele_freq_hist <- function(minor_allele_freq, col, main, adjust = 1.5, min_cov = 0){
	ks <- density(minor_allele_freq, adjust = adjust)
	peaks <- get_peaks(ks)
	print(peaks)
	hist(minor_allele_freq, col = col, freq = F, breaks = 20, xlab = 'minor allele coverage ratio', main = main, ylim = c(0, 12))
	lines(ks, lwd = 3)
	arrows(peaks[nrow(peaks), 'cov'], peaks[nrow(peaks), 'height'] + 1, peaks[nrow(peaks), 'cov'], peaks[nrow(peaks), 'height'] + 0.2, length = 0.1, lwd = 2)
	text(peaks[nrow(peaks), 'cov'], peaks[nrow(peaks), 'height'] + 1.4, round(peaks[nrow(peaks), 'cov'], 4), cex = 1)
}

figure_name = 'figures/het_autosomal_allele_supports/autosomal_coverage_ratio_Ocin2.png'
png(figure_name, units="in", width=5, height=5, res=300)
	plot_minor_allele_freq_hist(informative_A_snps$cov_minor / informative_A_snps$total_cov, 'grey', main = 'Ocin2', 2)
dev.off()
