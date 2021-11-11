
cov_tab <- read.table("tables/Afus_mean_coverage_table.tsv", header = T, check.names = F)

coverage_est_tab <- read.table("tables/resequencing_coverage_estimates.tsv", header = T, check.names = F)
row.names(coverage_est_tab) <- coverage_est_tab$ind

########################
# median covrage ratio
########################

males <- c('Afus1', 'BH3-2')
females <- coverage_est_tab$ind[!coverage_est_tab$ind %in% c('Afus1', 'BH3-2')]

cov_mat <- t(t(cov_tab[, coverage_est_tab$ind]) / coverage_est_tab[, 'diploid'])
cov_mat[is.na(cov_mat)] <- 0

male_covs <- rowMeans(cov_mat[, males])
female_covs <- apply(cov_mat[, females], 1, median, rm.na = T)

cov_tab$mf_cov_ratio <- log2(male_covs / female_covs)

# png('figures/Afus_comparing_male_and_population_asn.png')
#   plot( cov_tab$mf_cov_ratio[cov_tab$len > 10000] ~ log2(cov_tab$Afus1 / 50.52898)[cov_tab$len > 10000], xlim = c(-2, 2), ylim = c(-2, 2), xlab = 'log2 normalized ref male cov', ylab = 'log2 male / female cov ratio')
# dev.off()


scf_ks <- density(cov_tab$mf_cov_ratio, bw = "nrd0", adjust = 1, weights = cov_tab$len / sum(cov_tab$len))
second_deriv <- diff(sign(diff(scf_ks$y)))

peak_covs <- scf_ks$x[which(second_deriv == -2) + 1]
peak_heights <- scf_ks$y[which(second_deriv == -2) + 1]
# -0.65614821  0.02577876
# -0.372012 the local minima

cov_tab$asn_pop_means <- 'A'
cov_tab$asn_pop_means[cov_tab$mf_cov_ratio < -0.372012] <- 'X'

sum(cov_tab[cov_tab$asn_pop_means == 'A', 'len'])
sum(cov_tab[cov_tab$asn_pop_means == 'X', 'len'])

################

female_diffs <- rowSums(t(t(cov_tab[, females]) - coverage_est_tab[females, 'diploid']))
female_squares <- rowSums(t(t(cov_tab[, females]) - coverage_est_tab[females, 'diploid'])^2)
mean_norm_female_cov <- rowMeans(t(t(cov_tab[, females]) / coverage_est_tab[females, 'diploid']))

# male_squares_A <- rowSums(t(t(cov_tab[, males]) - coverage_est_tab[males, 'diploid'])^2)
# male_squares_X <- rowSums(t(t(cov_tab[, males]) - coverage_est_tab[males, 'haploid'])^2)

# cov_tab$A_squares <- female_squares + male_squares_A
# cov_tab$X_squares <- female_squares + male_squares_X
cov_tab$f_cov_var <- c(female_squares / 10)

# exploration of how it all behaves
# min_cov <- cov_tab$len > 3000
#
# hist(log10(female_squares[cov_tab$len > 20000]), breaks = 100)
#
# plot(log10(female_squares[min_cov]) ~ log10(cov_tab$len[min_cov]))
#
# subset <- sample(which(female_diffs < 2000 & min_cov), 2000)
# plot(female_diffs[subset] ~ log10(cov_tab$len[subset]), xlab = 'log10 scf len', ylab = 'sum of female diffs')
#
#
# plot(log10(male_squares_A[min_cov]) ~ log10(cov_tab$len[min_cov]))
# plot(log10(male_squares_X[min_cov]) ~ log10(cov_tab$len[min_cov]))
#
# plot(log10(male_squares_A[min_cov]) ~ log10(male_squares_X[min_cov]), xlim = c(-1, 5), ylim = c(-1, 5), xlab = 'log10 X squares', ylab = 'log10 A squares')
#
# plot(log10(male_squares_A[min_cov]) ~ log10(female_squares[min_cov]))
#
# plot(log10(male_squares_A[female_squares < 100]) ~ log10(male_squares_X[female_squares < 100]), xlim = c(-1, 5), ylim = c(-1, 5), xlab = 'log10 X squares', ylab = 'log10 A squares')
# plot(cov_tab$mf_cov_ratio[female_squares < 1000] ~ log10(female_squares[female_squares < 1000]), xlim = c(-1, 5), ylim = c(-1, 5))
#
# asn_stats <- function(var_filt, cov_filt, len_filt){
#   cov_tab$var_filt <- F
#   cov_tab$var_filt[cov_tab$f_cov_var < var_filt] <- T
#
#   cov_tab$cov_filt <- F
#   cov_tab$cov_filt[mean_norm_female_cov > cov_filt] <- T
#
#   cov_tab$len_filt <- F
#   cov_tab$len_filt[cov_tab$len > len_filt] <- T
#
#   assignable <- cov_tab$var_filt & cov_tab$cov_filt & cov_tab$len_filt
#
#   return(
#   c(sum(cov_tab[assignable, 'len']) / 1e6,
#     length(cov_tab[assignable, 'len']),
#     sum(cov_tab[assignable & cov_tab$asn_pop_means == 'A', 'len'], na.rm = T) / 1e6,
#     sum(cov_tab[assignable & cov_tab$asn_pop_means == 'X', 'len'], na.rm = T) / 1e6))
# }
#
# filt_asn_tab <- expand.grid(var = c(10, 100, 1000), cov = c(0, 0.3, 0.5), len = c(1000, 5000, 10000, 20000))
# filt_asn_tab[, c('assignable_sum', 'assignable_scf', 'A', 'X')] <- NA
#
# for (i in 1:nrow(filt_asn_tab)){
#   filt_asn_tab[i, c('assignable_sum', 'assignable_scf', 'A', 'X')] <- asn_stats(filt_asn_tab[i, 'var'], filt_asn_tab[i, 'cov'], filt_asn_tab[i, 'len'])
# }

# The coverage filter does nothing, is < 0.5 cov
cov_tab$var_filt <- F
cov_tab$var_filt[cov_tab$f_cov_var < 10] <- T

cov_tab$looser_var_filt <- F
cov_tab$looser_var_filt[cov_tab$f_cov_var < 100] <- T

cov_tab$len_filt <- F
cov_tab$len_filt[cov_tab$len > 10000] <- T

assignable <- cov_tab$var_filt & cov_tab$len_filt
c(sum(cov_tab[assignable, 'len']) / 1e6,
  length(cov_tab[assignable, 'len']),
  sum(cov_tab[assignable & cov_tab$asn_pop_means == 'A', 'len'], na.rm = T) / 1e6,
  sum(cov_tab[assignable & cov_tab$asn_pop_means == 'X', 'len'], na.rm = T) / 1e6)

cov_tab$asn_pop_means_filt <- NA

candidates <- cov_tab$len_filt & cov_tab$looser_var_filt
trusted <- cov_tab$len_filt & cov_tab$var_filt

cov_tab[candidates, 'asn_pop_means_filt'] <- paste0(cov_tab[candidates, 'asn_pop_means'], 'c')
cov_tab[trusted, 'asn_pop_means_filt'] <- cov_tab[trusted, 'asn_pop_means']

head(cov_tab[!is.na(cov_tab$asn_pop_means_filt),])

sapply(c("A", "Ac", "X", "Xc"), function(x) { sum(cov_tab[cov_tab$asn_pop_means_filt == x, 'len'], na.rm = T) } )

out_tab <- cov_tab[!is.na(cov_tab$asn_pop_means_filt), c('scf', 'len', 'mf_cov_ratio', 'f_cov_var', 'asn_pop_means_filt')]

out_tab[, c('mf_cov_ratio', 'f_cov_var')] <- round(out_tab[, c('mf_cov_ratio', 'f_cov_var')], 3)
colnames(out_tab)[5] <- 'chr'

write.table(out_tab, "tables/chrosmome_asn.tsv", row.names = F, col.names = T, quote = F, sep = '\t')
