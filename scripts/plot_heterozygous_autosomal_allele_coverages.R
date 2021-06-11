#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

ind <- args[1]
snp_tab <- read.table(args[2])
colnames(snp_tab) <- c('scf', 'pos', 'genotype', 'total_cov', 'ref_cov', 'alt_cov')

asn_tab <- read.table('tables/chrosmome_asn.tsv', header = T)
strlen <- max(nchar(asn_tab$scf))

names2tokens <- function(scf_name){
  sapply(strsplit(scf_name, '_'), function(x){ added_0s <- (strlen - (1 + sum(nchar(x)))); paste0(x[1], paste0(rep('0', added_0s), collapse = ''), x[2]) } )
}

rownames(asn_tab) <- names2tokens(asn_tab$scf)
snp_tab$chr <- asn_tab[names2tokens(snp_tab$scf), 'chr']

snp_tab$cov_minor <- apply(snp_tab[, c('ref_cov', 'alt_cov')], 1, min)

informative_A_snps <- snp_tab[snp_tab$chr == 'A' & snp_tab$genotype == '0/1', ]
informative_X_snps <- snp_tab[snp_tab$chr == 'X' & snp_tab$genotype == '0/1', ]

png(args[3])
  hist(informative_A_snps$cov_minor, col = 'blue', main = paste('Coverages supporting autosomal heterozygous alleles in', ind), xlab = 'Coverage')
  hist(informative_A_snps$total_cov - informative_A_snps$cov_minor, col = rgb(1, 0.05, 0.1, 0.75), add = T, breaks = 60)
  if (ind == 'BH3-2'){
    lines(c(11.3, 11.3), c(0, 1e6), lwd = 3, lty = 2)
    lines(c(18.5, 18.5), c(0, 1e6), lwd = 3, lty = 2)
  }

  legend('topright', pch = 20, col = c('blue', rgb(1, 0.05, 0.1, 0.75)), c('minor (paternal) allele', 'major (maternal) allele'))
dev.off()
