library('reldist')

########
# FUNCIONS
#######

get_peak <- function(ks, which_peak = 1){
  second_deriv <- diff(sign(diff(ks$y)))

  peak_covs <- ks$x[which(second_deriv == -2) + 1]
  peak_heights <- ks$y[which(second_deriv == -2) + 1]
  peak_covs[order(peak_heights, decreasing=T)][which_peak]
}

smoothing_ests <- function(cov_tab, ind, bw_technique, adjust){
  filt_quantile <- wtd.quantile(cov_tab[, ind], 0.98, weight = cov_tab[, 'len'])
  ind_tab <- cov_tab[cov_tab[, ind] < filt_quantile & cov_tab[, 'len'] > 20000, c('scf', 'len', ind)]
  ind_tab <- ind_tab[!is.na(ind_tab[, ind]), ]

  scf_ks <- density(ind_tab[, ind], bw = bw_technique, adjust = adjust, weights = ind_tab[, 'len'] / sum(ind_tab[, 'len']))
  est_1 <- get_peak(scf_ks, 1)
  est_2 <- get_peak(scf_ks, 2)
  fraction_of_sperm(est_1, est_2)
  c(est_2 - est_1, est_1, est_2, fraction_of_sperm(est_1, est_2), scf_ks$bw)
}

fraction_of_sperm <- function(X_cov, A_cov){ return( ((2 * X_cov) - A_cov) / X_cov  ) }

########
# MEANS
########
cov_means <- read.table('tables/Afus_mean_coverage_table.tsv', header = T, check.names = F)
# cov_medians <- read.table('tables/Afus_median_coverage_table.tsv', header = T, check.names = F)

smoothing_ests('BH3-2', 'SJ', 1)
# [1] 18.3265818 29.6746452  0.3807867
smoothing_ests('Afus1', 'SJ', 1)
# [1] 57.7197160 95.1939575  0.3507549

# testing stuff
# coverages <- ind_tab[, ind]
# scf_lengths <- ind_tab$len
# bw_technique <- 'SJ'
# adjust <- 1
#
# scf_ks <- density(coverages, bw = bw_technique, adjust = adjust, weights = scf_lengths / sum(scf_lengths))
# est_1 <- get_peak(scf_ks, 1)
# est_2 <- get_peak(scf_ks, 2)
# fraction_of_sperm(est_1, est_2)
# c(est_2 - est_1, est_1, est_2, fraction_of_sperm(est_1, est_2))

smoothing_ests(cov_means, 'BH3-2', 'SJ', 1)
# [1] 11.3033301 18.3541031 29.6574332  0.3841524  0.3058355
smoothing_ests(cov_medians, 'BH3-2', 'SJ', 1)
# [1] 10.97215804 18.01934702 28.99150505  0.39109014  0.04797985

smoothing_ests(cov_means, 'Afus1', 'SJ', 1)
# [1] 37.4459078 57.7044404 95.1503482  0.3510741  0.5062469
smoothing_ests(cov_medians, 'Afus1', 'SJ', 1)
# [1] 37.900072 57.015051 94.915122  0.335262  0.424832
