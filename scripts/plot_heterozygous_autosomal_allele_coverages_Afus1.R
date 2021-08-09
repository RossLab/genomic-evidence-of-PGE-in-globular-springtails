per_scf_cov_medians <- read.table('data/mapped_reads/Afus1_per_scf_cov_medians.tsv', col.names = c('scf', 'cov'))
asn_tab <- read.table('tables/chr_assignments_Afus1.tsv', header = T)

strlen <- max(nchar(c(asn_tab$scf, per_scf_cov_medians$scf)))

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(asn_tab) <- names2tokens(asn_tab$scf)
per_scf_cov_medians$chr <- asn_tab[names2tokens(per_scf_cov_medians$scf), 'chr']
per_scf_cov_medians$len <- asn_tab[names2tokens(per_scf_cov_medians$scf), 'len']

# per_scf_cov_medians <- per_scf_cov_medians[per_scf_cov_medians$cov_median < 60, ]
per_scf_cov_medians <- per_scf_cov_medians[!is.na(per_scf_cov_medians$chr), ]

autosomes <- per_scf_cov_medians$chr == 'A'
sex_chr <- per_scf_cov_medians$chr == 'X'

# plot(per_scf_cov_medians$len ~ per_scf_cov_medians$cov_median, pch = 20)
# points(per_scf_cov_medians$len[autosomes] ~ per_scf_cov_medians$cov_median[autosomes], col = 'yellow')
# points(per_scf_cov_medians$len[sex_chr] ~ per_scf_cov_medians$cov_median[sex_chr], col = 'green')

autosome_ks <- density(per_scf_cov_medians[autosomes, 'cov'], bw = "nrd0", adjust = 2.5, weights = per_scf_cov_medians$len[autosomes] / sum(per_scf_cov_medians$len[autosomes]))
second_deriv <- diff(sign(diff(autosome_ks$y)))

peak_covs <- autosome_ks$x[which(second_deriv == -2) + 1]
peak_heights <- autosome_ks$y[which(second_deriv == -2) + 1]

peak_covs
peak_heights
diploid = 95.03840


sex_ks <- density(per_scf_cov_medians[sex_chr, 'cov'], bw = "nrd0", adjust = 4, weights = per_scf_cov_medians$len[sex_chr] / sum(per_scf_cov_medians$len[sex_chr]))
plot(sex_ks)

second_deriv <- diff(sign(diff(autosome_ks$y)))

peak_covs <- sex_ks$x[which(second_deriv == -2) + 1]
peak_heights <- sex_ks$y[which(second_deriv == -2) + 1]

peak_covs
peak_heights
haploid = peak_covs[which.max(peak_heights)]

diploid - haploid

#!/usr/bin/env Rscript
ind <- 'Afus1'
snp_tab_A <- read.table('data/SNP_calls/freebayes_Afus_filt_sorted_Afus1_A.tsv')
snp_tab_X <- read.table('data/SNP_calls/freebayes_Afus_filt_sorted_Afus1_X.tsv')
colnames(snp_tab_A) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')
colnames(snp_tab_X) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')

snp_tab_A$cov_minor <- apply(snp_tab_A[, c('ref_cov', 'alt_cov')], 1, min)
snp_tab_X$cov_minor <- apply(snp_tab_X[, c('ref_cov', 'alt_cov')], 1, min)

ref_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/0', ]
  informative_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/1' & snp_tab_A$total_cov < 100, ]
informative_X_snps <- snp_tab_X[snp_tab_X$genotype == '0/0' & snp_tab_X$total_cov < 100, ]

pal <- c('black', rgb(0.8, 0.05, 0.1, 0.55), rgb(0.02, 0.45, 0.65, 0.55))

hist(informative_X_snps$total_cov, col = pal[1], main = paste('Coverages supporting autosomal heterozygous alleles in', ind), xlab = 'Coverage', breaks = 100, freq = F, xlim = c(0, 80), ylim = c(0, 0.05))
hist(informative_A_snps$cov_minor, col = pal[3], add = T, breaks = 30, freq = F, border = F)
hist(informative_A_snps$total_cov - informative_A_snps$cov_minor, col = pal[2], add = T, breaks = 30, freq = F, border = F)
if (ind == 'Afus1'){
  lines(c(37.6852, 37.6852), c(0, 1e6), lwd = 3, lty = 2)
  lines(c(57.3532, 57.3532), c(0, 1e6), lwd = 3, lty = 2)
}

legend('topright', pch = 20, col = pal, c('X chromosome allele', 'minor (paternal) allele', 'major (maternal) allele'), bty = 'n')
