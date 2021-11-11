library('weights')

get_peak <- function(ks, which_peak = 1){
  second_deriv <- diff(sign(diff(ks$y)))

  peak_covs <- ks$x[which(second_deriv == -2) + 1]
  peak_heights <- ks$y[which(second_deriv == -2) + 1]
  peak_covs[order(peak_heights, decreasing=T)][which_peak]
}

##################3
# JOIN PLOT MEANS #
##################3

cov_tab <- read.table("tables/Afus_mean_coverage_table.tsv", header = T, check.names = F)

coverage_est_tab <- data.frame(ind = c('Afus1', colnames(cov_tab)[4:15]), haploid = NA, diploid = NA)
row.names(coverage_est_tab) <- coverage_est_tab$ind

####
# COVERAGE SEXING PLOT
####
png('figures/Allacma_cov_estimates.png', width = 800, height = 600)
par(mfrow = c(3, 5))

for (ind in coverage_est_tab$ind){
  filt_quantile <- wtd.quantile(cov_tab[, ind], 0.98, weight = cov_tab[, 'len'])
  ind_tab <- cov_tab[cov_tab[, ind] < filt_quantile & cov_tab[, 'len'] > 20000, c('scf', 'len', ind)]
  ind_tab <- ind_tab[!is.na(ind_tab[, ind]), ]

  scf_ks <- density(ind_tab[, ind], bw = "SJ", adjust = 1, weights = ind_tab$len / sum(ind_tab$len))
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
  plot(scf_ks, main = ind, xlab = 'Coverage')
}

dev.off()

write.table(coverage_est_tab, "tables/resequencing_coverage_estimates.tsv", row.names = F, col.names = T, quote = F, sep = '\t')


####
# INDEIIDUAL COVERAGE PLOTS
####

prefix <- 'figures/mapping_coverages/'
suffix <- '_coverage_plot.png'

# Afus1
ind <- 'BH3-2'
figurename <- paste0(prefix, ind, suffix)

filt_quantile <- wtd.quantile(cov_tab[, ind], 0.98, weight = cov_tab[, 'len'])
ind_tab <- cov_tab[cov_tab[, ind] < filt_quantile & cov_tab[, 'len'] > 20000, c('scf', 'len', ind)]
ind_tab <- ind_tab[!is.na(ind_tab[, ind]), ]
#

adjust = 1
scf_ks <- density(ind_tab[, ind], bw = "SJ", adjust = adjust, weights = ind_tab$len / sum(ind_tab$len))
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

png(figurename, units="in", width=5, height=5, res=300)

  par(mar = c(4, 4, 1, 1) + 0.1)

  subset <- ind_tab[, ind] < 60
  wtd.hist(ind_tab[subset, ind], breaks = 60,
       freq = F, weight = ind_tab[subset, 'len'],
       col = 'grey', border = NA,
       main = '',
       xlim = c(0, 60), ylim = c(0, 0.30),
       xlab = 'Mean scaffold coverage')
	lines(scf_ks, lwd = 2)

  lines(c(coverage_est_tab[ind, 'haploid'], coverage_est_tab[ind, 'haploid']), c(0, 1000), lty = 2, lwd = 2)
  lines(c(coverage_est_tab[ind, 'diploid'], coverage_est_tab[ind, 'diploid']), c(0, 1000), lty = 2, lwd = 2)
  # mtext(paste0(round(coverage_est_tab[ind, 'haploid'], 1),'x'), 3, padj = 2, at = coverage_est_tab[ind, 'haploid'] + 5, cex = 1)
  # mtext(paste0(round(coverage_est_tab[ind, 'diploid'], 1),'x'), 3, padj = 2, at = coverage_est_tab[ind, 'diploid'] + 5, cex = 1)


dev.off()

ind = 'Afus1'
figurename <- paste0(prefix, ind, suffix)

filt_quantile <- wtd.quantile(cov_tab[, ind], 0.98, weight = cov_tab[, 'len'])
ind_tab <- cov_tab[cov_tab[, ind] < filt_quantile & cov_tab[, 'len'] > 20000, c('scf', 'len', ind)]

adjust = 1
scf_ks <- density(ind_tab[, ind], bw = "SJ", adjust = adjust, weights = ind_tab$len / sum(ind_tab$len))

png(figurename, units="in", width=5, height=4, res=300)

  par(mar = c(4, 4, 1, 1) + 0.1)

  subset <- ind_tab[, ind] < filt_quantile
  wtd.hist(ind_tab[subset, ind], breaks = 60,
       freq = F, weight = ind_tab[subset, 'len'],
       col = 'grey', border = NA,
       main = '',
       xlim = c(0, 140), ylim = c(0, 0.20),
       xlab = 'Mean scaffold coverage')
	lines(scf_ks, lwd = 2)

  lines(c(coverage_est_tab[ind, 'haploid'], coverage_est_tab[ind, 'haploid']), c(0, 1000), lty = 2, lwd = 2)
  lines(c(coverage_est_tab[ind, 'diploid'], coverage_est_tab[ind, 'diploid']), c(0, 1000), lty = 2, lwd = 2)

dev.off()


# ###########
# # MEDIANS #
# ###########
#
# median_cov_tab <- read.table("tables/Afus_median_coverage_table.tsv", header = T, check.names = F)
# strlen <- max(nchar(c(cov_tab$scf, median_cov_tab$scf)))
#
# names2tokens <- function(scf_name){
#   sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
# }
#
# rownames(median_cov_tab) <- names2tokens(median_cov_tab$scf)
# median_cov_tab_reduced <- median_cov_tab[names2tokens(cov_tab$scf), ]
#
# cov_tab <- median_cov_tab_reduced


#########
# Ocin2 #
#########

# this table also contain Ocin2 coverages in "male_coverage" column
asn_tab <- read.table('tables/chr_assignments_Ocin1.tsv', header = T)

ind = 'Ocin2'
figurename <- paste0(prefix, ind, suffix)

filt_quantile <- wtd.quantile(asn_tab[, 'male_coverage'], 0.98, weight = asn_tab[, 'len'])
ind_tab <- asn_tab[asn_tab[, 'male_coverage'] < filt_quantile & asn_tab[, 'len'] > 20000, c('scf', 'len', 'male_coverage')]

adjust = 1
scf_ks <- density(ind_tab[, 'male_coverage'], bw = "nrd0", adjust = adjust, weights = ind_tab$len / sum(ind_tab$len))
cov_2n <- get_peak(scf_ks)
# [1] 101.1228
cov_1n <- get_peak(scf_ks, 2)
# [1] 53.4376

png(figurename, units="in", width=5, height=4, res=300)

  par(mar = c(4, 4, 1, 1) + 0.1)

  subset <- ind_tab[, 'male_coverage'] < filt_quantile
  wtd.hist(ind_tab[subset, 'male_coverage'], breaks = 60,
       freq = F, weight = ind_tab[subset, 'len'],
       col = 'grey', border = NA,
       main = '',
       xlim = c(0, 150), ylim = c(0, 0.05),
       xlab = 'Mean scaffold coverage')
  lines(scf_ks, lwd = 2)


  lines(c(cov_1n, cov_1n), c(0, 1000), lty = 2, lwd = 2)
  lines(c(cov_2n, cov_2n), c(0, 1000), lty = 2, lwd = 2)

dev.off()
