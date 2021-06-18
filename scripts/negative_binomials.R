allele_1 <- rnbinom(200000, size = 15, mu = 14)
allele_2 <- rnbinom(200000, size = 15, mu = 14)

total <- allele_1 +  allele_2

minor_cov <- apply(data.frame(cov1 = allele_1, cov2 = allele_2), 1, min)
pal <- c('black', rgb(0.02, 0.45, 0.65, 0.55), rgb(0.8, 0.05, 0.1, 0.55))
cex_legend = 1.5

png('figures/simulated_allele_cov_female.png')

hist(allele_1, col = pal[1], main = '', xlab = 'Coverage', breaks = 50, freq = F, cex.axis = cex_legend, cex.lab = cex_legend, ylim = c(0, 0.1))
hist(minor_cov, col = pal[2], add = T, breaks = 30, freq = F, border = F)
hist(total - minor_cov, col = pal[3], add = T, breaks = 60, freq = F, border = F)

legend('topright', pch = 20, col = pal, c('real allele coverages', 'minor alleles', 'major alleles'), bty = 'n', cex = cex_legend, title = 'simulated coverages')

dev.off()

allele_1 <- rnbinom(200000, size = 15, mu = 10)
allele_2 <- rnbinom(200000, size = 15, mu = 19)

total <- allele_1 +  allele_2
minor_cov <- apply(data.frame(cov1 = allele_1, cov2 = allele_2), 1, min)

png('figures/simulated_allele_cov_male.png')

hist(allele_1, col = pal[1], main = ' ', xlab = 'Coverage', breaks = 50, freq = F, cex.axis = cex_legend, cex.lab = cex_legend, ylim = c(0, 0.1))
hist(minor_cov, col = pal[2], add = T, breaks = 30, freq = F, border = F)
hist(total - minor_cov, col = pal[3], add = T, breaks = 60, freq = F, border = F)

legend('topright', pch = 20, col = pal, c('real paternal coverages', 'minor alleles', 'major alleles'), bty = 'n', cex = cex_legend, title = 'simulated coverages')


dev.off()
