# 'data/mapped_reads/per_scf_cov_medians_Afus1.tsv'
BH3_medians <- read.table('data/mapped_reads/per_scf_cov_medians_BH3-2.tsv', col.names = c('scf', 'cov_median'))

asn_tab <- read.table('tables/chr_assignments_Afus1.tsv', header = T)
strlen <- max(nchar(c(asn_tab$scf, BH3_medians$scf)))

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(asn_tab) <- names2tokens(asn_tab$scf)
BH3_medians$chr <- asn_tab[names2tokens(BH3_medians$scf), 'chr']
BH3_medians$len <- asn_tab[names2tokens(BH3_medians$scf), 'len']

BH3_medians <- BH3_medians[BH3_medians$cov_median < 60, ]
BH3_assigned_medians <- BH3_medians[!is.na(BH3_medians$chr), ]

autosomes <- BH3_assigned_medians$chr == 'A'
sex_chr <- BH3_assigned_medians$chr == 'X'

# plot(BH3_assigned_medians$len ~ BH3_assigned_medians$cov_median, pch = 20)
# points(BH3_assigned_medians$len[autosomes] ~ BH3_assigned_medians$cov_median[autosomes], col = 'yellow')
# points(BH3_assigned_medians$len[sex_chr] ~ BH3_assigned_medians$cov_median[sex_chr], col = 'green')

BH3_autosome_ks <- density(BH3_assigned_medians[autosomes, 'cov_median'], bw = "nrd0", adjust = 2.5, weights = BH3_assigned_medians$len[autosomes] / sum(BH3_assigned_medians$len[autosomes]))
second_deriv <- diff(sign(diff(BH3_autosome_ks$y)))

peak_covs <- BH3_autosome_ks$x[which(second_deriv == -2) + 1]
peak_heights <- BH3_autosome_ks$y[which(second_deriv == -2) + 1]

peak_covs
peak_heights
diploid = 29.49488


BH3_sex_ks <- density(BH3_assigned_medians[sex_chr, 'cov_median'], bw = "nrd0", adjust = 3, weights = BH3_assigned_medians$len[sex_chr] / sum(BH3_assigned_medians$len[sex_chr]))
second_deriv <- diff(sign(diff(BH3_autosome_ks$y)))

peak_covs <- BH3_sex_ks$x[which(second_deriv == -2) + 1]
peak_heights <- BH3_sex_ks$y[which(second_deriv == -2) + 1]

peak_covs
peak_heights
haploid = 19.467778

diploid - haploid
