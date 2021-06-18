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

pal <- c('black', rgb(0.02, 0.45, 0.65, 0.55), rgb(0.8, 0.05, 0.1, 0.55))
cex_legend <- 1.4
# xlim <- c(0, 50)
# ylim <- c(0, 0.1)

png(args[3])
  hist(informative_A_snps$cov_minor, col = pal[2], main = paste('Coverages supporting alleles in', ind), xlab = 'Coverage', freq = F, border = F, cex.axis = cex_legend, cex.lab = cex_legend)
  hist(informative_A_snps$total_cov - informative_A_snps$cov_minor, col = pal[3], add = T, breaks = 60, freq = F, border = F)
  if (ind == 'BH3-2'){
    lines(c(10.03, 10.03), c(0, 1e6), lwd = 3, lty = 2)
    lines(c(19.47, 19.47), c(0, 1e6), lwd = 3, lty = 2)
  }

  legend('topright', pch = 20, col = pal[c(2, 3)], c('minor allele', 'major allele'))
dev.off()
