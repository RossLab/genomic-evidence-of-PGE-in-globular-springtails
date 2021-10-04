allele_1 <- rnbinom(200000, size = 30, mu = 52)
allele_2 <- rnbinom(200000, size = 30, mu = 52)

total <- allele_1 +  allele_2

minor_cov <- apply(data.frame(cov1 = allele_1, cov2 = allele_2), 1, min)

paste("Fraction of assigned: ", mean(minor_cov == allele_1))
source('scripts/load_palette.R')
cex_legend = 1 #0.95

source('scripts/fixed_bin_historgram.R')

coverage_data <- list(maternal = allele_1, minor = minor_cov, major = total - minor_cov)

figure_name <- 'figures/simulated_allele_cov_female.png'
png(figure_name, units="in", width=5, height=5, res=300)

par(mar = c(4, 4, 1, 1) + 0.1)
fixed_bin_histogram(coverage_data, pal, xlab = 'Coverage support', bins = 50, freq = F, xlim = c(0, 120), default_legend = F)
# legend('topright', pch = 20, col = pal, c('real coverage distribution', 'minor alleles', 'major alleles'), bty = 'n', cex = cex_legend, title = 'simulated coverages')
legend('topright', pch = 20, col = pal, c('X-linked homozygous', 'minor autosomal', 'major autosomal'), bty = 'n', cex = cex_legend, title = 'allele coverage supports')

dev.off()

allele_1 <- rnbinom(200000, size = 15, mu = 11.3)
allele_2 <- rnbinom(200000, size = 15, mu = 18.35)

total <- allele_1 +  allele_2
minor_cov <- apply(data.frame(cov1 = allele_1, cov2 = allele_2), 1, min)
coverage_data <- list(maternal = allele_2, minor = minor_cov, major = total - minor_cov)

prop.test(sum(minor_cov != allele_1), 200000, conf.level=0.95, correct = FALSE)
# paste("Fraction of assigned: ",)

figure_name <- 'figures/simulated_allele_cov_male.png'
png(figure_name, units="in", width=5, height=5, res=300)

par(mar = c(4, 4, 1, 1) + 0.1)
fixed_bin_histogram(coverage_data, pal, xlab = 'Coverage support', bins = 50, freq = F, xlim = c(0, 50), default_legend = F)
# legend('topright', pch = 20, col = pal, c('real maternal coverages', 'minor alleles', 'major alleles'), bty = 'n', cex = cex_legend, title = 'simulated coverages')

dev.off()
