barplot <- function(bar_heights, bar_positions, width = 0.5, xlim = NA, ylim = NA, font_size = 1, xlab = '', ylab = '', add = F, main = '', col = 'deepskyblue', yaxis = T){
	if ( any(is.na(xlim)) ){
		xlim <- range(bar_positions) + c(-0.8, 0.8)
	}
	if ( any(is.na(ylim)) ){
		ylim <- range(bar_heights)
	}
	if ( add == F ){
		plot(bar_heights, type="n", xlab=xlab, ylab=ylab,
			ylim=ylim, xlim=xlim, xaxt="n", yaxt = 'n', bty = 'n',
			cex.lab=font_size, cex.main=font_size, cex.sub=font_size, main = main)
		if ( yaxis ){
			axis(2, cex.axis=font_size)
		}
	}
	for ( i in 1:length(bar_heights)){
		rect(bar_positions[i] - width, 0, bar_positions[i] + width, bar_heights[i], col = col, border = F)
	}
}

all_sims <- read.table('tables/power_analysis_complete.tsv', header = T, sep = '\t')
all_sims <- all_sims[order(all_sims$theta),]

pdf('figures/power_analysis_proportional.pdf')
	par(mfrow = c(4, 3))
	for(X in unique(all_sims$X)){
		is_X <- all_sims$X == X
		for(mat in unique(all_sims$maternal)){
			is_coverage <- all_sims$maternal == mat
			# for(theta in unique(all_sims$theta)){
				# is_theta <- all_sims$theta == theta
				ylim <- c(-0.05, 0.55)
				to_plot <- is_X & is_coverage # & is_theta
				sim_sperm <- all_sims$simulated_sperm[to_plot]
				kmer_sperm <- all_sims$two_tissue_sperm[to_plot]
				mapping_sperm <- all_sims$mapping_sperm[to_plot]

				if ( mat == 10 ){
					# bottom, left, top and right margins
					par(mar=c(0.2,3.1,0.1,0.6))
					barplot(sim_sperm, 1:sum(to_plot), ylim = ylim, yaxis = T)
				} else {
					par(mar=c(0.2,0,0.1,0.2))
					barplot(sim_sperm, 1:sum(to_plot), ylim = ylim, yaxis = F)
				}
				barplot(kmer_sperm, 1:sum(to_plot) - 0.25, add = T, width = 0.25, col = 'deeppink3')
				barplot(mapping_sperm, 1:sum(to_plot) + 0.25, add = T, width = 0.25, col = 'darkgoldenrod1')
			# }
		}
	}
dev.off()

# this image was then labeled and annotated in InkScape
all_sims$true_positive <- all_sims$CI_l_prop > 0.5 & all_sims$simulated_sperm > 0

print('True positives: ')
sum(all_sims$true_positive)
print('Decomposed by simulated parameters: ')
print('Fraction of sperm: ')
table(all_sims[all_sims$true_positive, 'simulated_sperm']) / 48
print('Number of X chrmosomes (out of 20): ')
table(all_sims[all_sims$true_positive, 'X']) / 60
print('heterozygosity (theta): ')
table(all_sims[all_sims$true_positive, 'theta']) / 60
print('coverage (1n maternal): ')
table(all_sims[all_sims$true_positive, 'maternal']) / 80
