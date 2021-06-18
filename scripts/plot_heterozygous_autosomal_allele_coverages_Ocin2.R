coverages <- read.table('data/mapped_reads/per_scf_cov_medians_Ocin2.tsv', col.names = c('scf', 'cov'))
# hist(coverages$cov[coverages$cov < 150], breaks = 150)
# 30 - 130 is reasonable cov range

ind <- 'Ocin2'
snp_tab_A <- read.table('data/SNP_calls/freebayes_Ocin2_filt_sorted_A_Ocin2.tsv')
snp_tab_X <- read.table('data/SNP_calls/freebayes_Ocin2_filt_sorted_X_Ocin2.tsv')
colnames(snp_tab_A) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')
colnames(snp_tab_X) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')


snp_tab_A$cov_minor <- apply(snp_tab_A[, c('ref_cov', 'alt_cov')], 1, min)
snp_tab_X$cov_minor <- apply(snp_tab_X[, c('ref_cov', 'alt_cov')], 1, min)

informative_A_snps <- snp_tab_A[snp_tab_A$genotype == '0/1' & snp_tab_A$total_cov < 120 & snp_tab_A$total_cov > 25, ]
informative_X_snps <- snp_tab_X[snp_tab_X$genotype == '1/1' & snp_tab_X$total_cov < 120 & snp_tab_X$total_cov > 120, ]

pal <- c('black', rgb(0.8, 0.05, 0.1, 0.55), rgb(0.02, 0.45, 0.65, 0.55))

hist( - informative_X_snps$cov_minor, col = pal[1], main = paste('Coverages supporting autosomal heterozygous alleles in', ind), xlab = 'Coverage', breaks = 60, freq = F, xlim = c(0, 120), ylim = c(0, 0.05))
hist(informative_A_snps$cov_minor, col = pal[2], add = T, breaks = 30, freq = F, border = F)
hist(informative_A_snps$total_cov - informative_A_snps$cov_minor, col = pal[3], add = T, breaks = 60, freq = F, border = F)

legend('topright', pch = 20, col = pal, c('X chromosome allele', 'minor allele', 'major allele'), bty = 'n')

# hist(informative_X_snps)

hist(informative_A_snps$cov_minor / informative_A_snps$total_cov)
