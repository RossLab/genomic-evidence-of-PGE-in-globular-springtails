
cluster_string <- 'qsub -o logs -e logs -cwd -N power_analysis -V -pe smp64 4 -b yes -l h="bigbang"'


heterozygosity <- c(0.0001, 0.001, 0.003, 0.005)
X_chromosomes <- c(1, 2, 5, 10)
fraction_of_sperm <- c(0.0, 0.01, 0.05, 0.1, 0.25, 0.50)
sequencing_depth <- c(10, 15, 25)

parameter_combinations <- expand.grid(heterozygosity = heterozygosity, X_chromosomes = X_chromosomes, fraction_of_sperm = fraction_of_sperm, sequencing_depth = sequencing_depth)

parameter_combinations$paternal_depth <- parameter_combinations$sequencing_depth * (1 - parameter_combinations$fraction_of_sperm )

commands <-	paste0(cluster_string,
	                 ' "scripts/run_power_analysis_replicate.sh ',
					              format(parameter_combinations$heterozygosity, scientific = F), ' ',
						  		      parameter_combinations$X_chromosomes, ' ',
							  			  parameter_combinations$sequencing_depth, ' ',
								  		  parameter_combinations$paternal_depth, '"')

writeLines(commands, 'scripts/power_analysis_commnads.sh')
