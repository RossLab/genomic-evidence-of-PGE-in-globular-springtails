barplot <- function(bar_heights, bar_positions, width = 0.5, xlim = NA, ylim = NA, font_size = 1, xlab = '', ylab = '', add = F, main = '', col = 'deepskyblue'){
	if ( any(is.na(xlim)) ){
		xlim <- range(bar_positions)
	}
	if ( any(is.na(ylim)) ){
		ylim <- range(bar_heights)
	}
	if ( add == F ){
		plot(bar_heights, type="n", xlab=xlab, ylab=ylab,
		ylim=ylim, xlim=xlim, xaxt="n", bty = 'n',
		cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size, main = main)
	}
	for ( i in 1:length(bar_heights)){
		rect(bar_positions[i] - width, 0, bar_positions[i] + width, bar_heights[i], col = col, border = F)
	}
}

all_sims <- read.table('tables/power_analysis_complete.tsv', header = T, sep = '\t')

pdf('figures/power_analysis.pdf')
	par(mfrow = c(4, 3))
	par(mar=c(0.6,3.6,1.1,0.6))
	for(X in unique(all_sims$X)){
		is_X <- all_sims$X == X
		for(mat in unique(all_sims$maternal)){
			is_coverage <- all_sims$maternal == mat
			ylim <- c(-0.1, 0.7)
			barplot(all_sims$simulated_sperm[is_X & is_coverage], 1:sum(is_X & is_coverage), ylim = ylim)
			barplot(all_sims$two_tissue_sperm[is_X & is_coverage], 1:sum(is_X & is_coverage) - 0.25, add = T, width = 0.25, col = 'deeppink3')
			barplot(all_sims$mapping_sperm[is_X & is_coverage], 1:sum(is_X & is_coverage) + 0.25, add = T, width = 0.25, col = 'darkgoldenrod1')
		}
	}
dev.off()

# this image was then labeled and annotated in InkScape
