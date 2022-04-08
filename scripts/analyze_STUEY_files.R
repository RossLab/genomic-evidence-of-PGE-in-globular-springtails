# I removed all quotemarks beforehand with sed -i 's/\"//g'

get_scaffold_names <- function(file){
	sapply(file, function(x){ x[1] } )
}

get_median_coverages <- function(file){
	as.numeric(sapply(file, function(x){ x[2] } ))
}

get_site_list <- function(file){
	lapply(file, function(line){matrix(as.numeric(line[c(-1, -2)]), ncol = 3, byrow = T)} )
}

# USAGE
# file <- strsplit(readLines('data/STUEYdata/BH3n2_Xc_siteCOVpIN.tsv'), '\t')
#
# scaffolds <- get_scaffold_names(file)
# median_covs <- get_median_coverages(file)
# all_sites <-  get_site_list(file)
#
# sites_per_sfv <- sapply(all_sites, function(scf_tab){ nrow(scf_tab) } ) # get number of sites per scf, to be able to reconstruct back which comes from where
# all_sites <- do.call(rbind, all_sites) # convert to matrix
#
# hist(all_sites[all_sites[, 3] == 1, 2], breaks = 20)

get_all_sites <- function(filename){
	file <- strsplit(readLines(filename), '\t')
	scaffolds <- get_scaffold_names(file)
	# median_covs <- get_median_coverages(file)
	all_sites <-  get_site_list(file)

	sites_per_sfv <- sapply(all_sites, function(scf_tab){ nrow(scf_tab) } ) # get number of sites per scf, to be able to reconstruct back which comes from where
	all_sites <- as.data.frame(do.call(rbind, all_sites)) # convert to matrix
	colnames(all_sites) <- c('median', 'cov_p', 'filter')
	all_sites$scf <- rep(scaffolds, times = sites_per_sfv)
	return(all_sites)
}

source('scripts/get_peaks.R')

pal <- c('#00BFC4', '#F8766D')

BH3_2_A <- get_all_sites('data/STUEYdata/BH3n2_A_siteCOVpIN.tsv')
BH3_2_X <- get_all_sites('data/STUEYdata/BH3n2_X_siteCOVpIN.tsv')

Afus1_2_A <- get_all_sites('data/STUEYdata/Afus1_A_siteCOVpIN.tsv')
Afus1_2_X <- get_all_sites('data/STUEYdata/Afus1_X_siteCOVpIN.tsv')

plot_minor_allele_freq_hist <- function(tab, col, main, adjust = 1.5, min_cov = 0){
	minor_allele_freq <- tab[(tab$median * tab$cov_p) > min_cov & tab$filter == 0, 'cov_p']
	ks <- density(minor_allele_freq, adjust = adjust)
	peaks <- get_peaks(ks)
	print(peaks)
	hist(minor_allele_freq, col = col, freq = F, breaks = 20, xlab = 'minor allele frequencies', main = main)
	lines(ks, lwd = 3)
	arrows(peaks[1, 'cov'], peaks[1, 'height'] + 2, peaks[1, 'cov'], peaks[1, 'height'] + 0.5, length = 0.1, lwd = 3)
	text(peaks[1, 'cov'], peaks[1, 'height'] + 2.3, round(peaks[1, 'cov'], 4), cex = 1.5)
}

p2f_h <- function(p){
	(1 - (2 * p)) / (1 - p)
}

plot_minor_allele_freq_hist(BH3_2_A, pal[2], 'BH3-2 Autosomes')
plot_minor_allele_freq_hist(BH3_2_X, pal[1], 'BH3-2 X-linked scfs')

plot_minor_allele_freq_hist(Afus1_2_A, pal[2], 'Afus1 Autosomes', 2, 1)
plot_minor_allele_freq_hist(Afus1_2_X, pal[1], 'Afus1 X-linked scfs', 2, 1)

tab = BH3_2_X
min_cov = 0
minor_allele_freq <- tab[(tab$median * tab$cov_p) > min_cov & tab$filter == 0, 'cov_p']
hist_for_logplot <- hist(minor_allele_freq, col = col, freq = F, breaks = 20, xlab = 'minor allele frequencies', main = main)
hist_for_logplot$counts <- log10(hist_for_logplot$counts)
plot(hist_for_logplot, freq = T, col = col)

hist(tab$median * tab$cov_p)
