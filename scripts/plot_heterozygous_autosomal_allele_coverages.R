#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

ind <- args[1]
snp_tab <- read.table(args[2])
figure_name = args[3]

colnames(snp_tab) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')

asn_tab <- read.table('tables/chr_assignments_Afus1.tsv', header = T)
strlen <- max(nchar(asn_tab$scf))

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(asn_tab) <- names2tokens(asn_tab$scf)
snp_tab$chr <- asn_tab[names2tokens(snp_tab$scf), 'chr']

snp_tab$cov_minor <- apply(snp_tab[, c('ref_cov', 'alt_cov')], 1, min)

informative_A_snps <- snp_tab[snp_tab$chr == 'A' & snp_tab$genotype == '0/1', ]
informative_X_snps <- snp_tab[snp_tab$chr == 'X' & snp_tab$genotype == '0/1', ]

# X = informative_X_snps$total_cov,
coverage_data <- list(A_minor = informative_A_snps$cov_minor, A_major = informative_A_snps$total_cov - informative_A_snps$cov_minor)
pal <- c(rgb(0.8, 0.05, 0.1, 0.55), rgb(0.02, 0.45, 0.65, 0.55)) #rgb(0, 0, 0, 0.7),
# source('scripts/load_palette.R')

main <- ''# paste('Coverages supporting autosomal heterozygous alleles in', ind)
cex_legend <- 0.75
xlim <- c(0, 40)
# ylim <- c(0, 0.1)

source('scripts/fixed_bin_historgram.R')

# png(figure_name)
png(figure_name, units="in", width=5, height=5, res=300)

  # 'c(bottom, left, top, right)'
  par(mar = c(4, 4, 1, 1) + 0.1)
  fixed_bin_histogram(coverage_data, pal, main = main, xlab = 'Coverage support', bins = 50, freq = F, default_legend = F, xlim = xlim) # ylim = ylim,
  legend('topright', pch = 20, col = pal, c('minor allele', 'major allele'), bty = 'n', cex = cex_legend)
dev.off()
