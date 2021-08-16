library('reldist')

#########
# MEANS #
#########

cov_tab <- read.table("tables/Afus_mean_coverage_table.tsv", header = T, check.names = F)

coverage_est_tab <- data.frame(ind = c('Afus1', colnames(cov_tab)[4:15]), haploid = NA, diploid = NA)
row.names(coverage_est_tab) <- coverage_est_tab$ind

####
# COVERAGE SEXING PLOT
####
pdf('figures/Allacma_cov_estimates.pdf', width = 8, height = 6)
par(mfrow = c(3, 5))

for (ind in coverage_est_tab$ind){
  filt_quantile <- wtd.quantile(cov_tab[, ind], 0.95, weight = cov_tab[, 'len'])
  ind_tab <- cov_tab[cov_tab[, ind] < filt_quantile & cov_tab[, 'len'] > 20000, c('scf', 'len', ind)]
  ind_tab <- ind_tab[!is.na(ind_tab[, ind]), ]

  if ( ind %in% c('Afus1', 'BH3-2')){
    adjust = 1.5
  } else {
    adjust = 3
  }

  scf_ks <- density(ind_tab[, ind], bw = "nrd0", adjust = adjust, weights = ind_tab$len / sum(ind_tab$len))
  second_deriv <- diff(sign(diff(scf_ks$y)))

  peak_covs <- scf_ks$x[which(second_deriv == -2) + 1]
  peak_heights <- scf_ks$y[which(second_deriv == -2) + 1]

  if ( ind %in% c('Afus1', 'BH3-2')){
    # male
    est_1 <- peak_covs[which.max(peak_heights)]
    est_2 <- peak_covs[peak_heights!=max(peak_heights)][which.max( peak_heights[peak_heights!=max(peak_heights)] )]
    if (est_1 > est_2){
        coverage_est_tab[ind, 'haploid'] <- est_2
        coverage_est_tab[ind, 'diploid'] <- est_1
    } else {
        coverage_est_tab[ind, 'haploid'] <- est_1
        coverage_est_tab[ind, 'diploid'] <- est_2
    }
  } else {
    # We expect only one big peak in females
    coverage_est_tab[ind, 'diploid'] <- peak_covs[which.max(peak_heights)]
  }

  # plot for inspecting
  plot(scf_ks, main = ind)
}

dev.off()

# write.table(coverage_est_tab, "tables/resequencing_coverage_estimates.tsv", row.names = F, col.names = T, quote = F, sep = '\t')


####
# INDEIIDUAL COVERAGE PLOTS
####

prefix <- 'figures/mapping_coverages/'
suffix <- '_coverage_plot.png'

# Afus1
ind <- 'BH3-2'
figurename <- paste0(prefix, ind, suffix)

get_peak <- function(ks){
  second_deriv <- diff(sign(diff(ks$y)))

  peak_covs <- ks$x[which(second_deriv == -2) + 1]
  peak_heights <- ks$y[which(second_deriv == -2) + 1]
  peak_covs[which.max(peak_heights)]
}

png(figurename)

  filt_quantile <- wtd.quantile(cov_tab[, ind], 0.95, weight = cov_tab[, 'len'])
  ind_tab <- cov_tab[cov_tab[, ind] < filt_quantile & cov_tab[, 'len'] > 20000, c('scf', 'len', ind)]
  ind_tab <- ind_tab[!is.na(ind_tab[, ind]), ]
#

  adjust = 1
  scf_ks <- density(ind_tab[, ind], bw = "nrd0", adjust = adjust, weights = ind_tab$len / sum(ind_tab$len))
  second_deriv <- diff(sign(diff(scf_ks$y)))

  peak_covs <- scf_ks$x[which(second_deriv == -2) + 1]
  peak_heights <- scf_ks$y[which(second_deriv == -2) + 1]

  est_1 <- peak_covs[which.max(peak_heights)]
  est_2 <- peak_covs[peak_heights!=max(peak_heights)][which.max( peak_heights[peak_heights!=max(peak_heights)] )]
  if (est_1 > est_2){
      coverage_est_tab[ind, 'haploid'] <- est_2
      coverage_est_tab[ind, 'diploid'] <- est_1
  } else {
      coverage_est_tab[ind, 'haploid'] <- est_1
      coverage_est_tab[ind, 'diploid'] <- est_2
  }

	plot(scf_ks, main = ind)

  if ( ind %in% c('Afus1', 'BH3-2')){
    lines(c(coverage_est_tab[ind, 'haploid'], coverage_est_tab[ind, 'haploid']), c(0, 1000), lty = 2)
    lines(c(coverage_est_tab[ind, 'diploid'], coverage_est_tab[ind, 'diploid']), c(0, 1000), lty = 2)
    text(0, 0.04, round(coverage_est_tab[ind, 'diploid'] - coverage_est_tab[ind, 'haploid'], 2))
    text(0, 0.05, round(coverage_est_tab[ind, 'haploid'], 2))
  }


dev.off()


###########
# MEDIANS #
###########

median_cov_tab <- read.table("tables/Afus_median_coverage_table.tsv", header = T, check.names = F)
strlen <- max(nchar(c(cov_tab$scf, median_cov_tab$scf)))

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(median_cov_tab) <- names2tokens(median_cov_tab$scf)
median_cov_tab_reduced <- median_cov_tab[names2tokens(cov_tab$scf), ]

cov_tab <- median_cov_tab_reduced
