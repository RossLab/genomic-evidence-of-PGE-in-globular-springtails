two_tissue_model_files <- list.files(path = "data/simulations", pattern = 'two_tissue_model_summary_tsv', recursive = T, full.names = T)


all_sims <- as.data.frame(t(sapply(two_tissue_model_files, read.table, header = T)))
sim_directories <- as.vector(sapply(row.names(all_sims), function(x){ strsplit(x, split = '/')[[1]][3] }))
sim_info <- strsplit(sim_directories, "_")

all_sims <- as.data.frame(sapply(all_sims, as.numeric))
all_sims$sim_id <- sapply(sim_info, function(x){ as.numeric(x[2]) })
all_sims$theta <- sapply(sim_info, function(x){ as.numeric(sub('.', '', x[3])) })
all_sims$X <- sapply(sim_info, function(x){ as.numeric(sub('.', '', x[4])) })
all_sims$maternal <- sapply(sim_info, function(x){ as.numeric(sub('.', '', x[5])) })
all_sims$paternal <- sapply(sim_info, function(x){ as.numeric(sub('.', '', x[6])) })

row.names(all_sims) <- all_sims$sim_id

mapping_kernel_smoothing_files <- list.files(path = "data/simulations", pattern = 'mapping_coverage_summary_tsv', recursive = T, full.names = T)
mapping_summary_table <- as.data.frame(t(sapply(mapping_kernel_smoothing_files, read.table, header = T)))

sim_directories <- as.vector(sapply(row.names(mapping_summary_table), function(x){ strsplit(x, split = '/')[[1]][3] }))
sim_info <- strsplit(sim_directories, "_")

mapping_summary_table <- as.data.frame(sapply(mapping_summary_table, as.numeric))
mapping_summary_table$sim_id <- sapply(sim_info, function(x){ as.numeric(x[2]) })

all_sims <- merge(all_sims, mapping_summary_table)

all_sims$simulated_sperm <- (all_sims$maternal - all_sims$paternal) / all_sims$maternal

two_tissue_cov_differences <- sapply(1:nrow(all_sims), function(i){ abs(all_sims$kmercov1[i] - all_sims$kmercov2[i]) })
two_tissue_1n_cov <- sapply(1:nrow(all_sims), function(i){ min(c(all_sims$kmercov1[i], all_sims$kmercov2[i])) })
# all_sims$two_tissue_sperm <- 1 - (two_tissue_cov_differences / two_tissue_1n_cov)
all_sims$two_tissue_sperm <- 1 - ((all_sims$kmercov2 - all_sims$kmercov1) / all_sims$kmercov1)

all_sims$mapping_sperm <-  1 - ((all_sims$mapping_disomic - all_sims$mapping_monosomic) / all_sims$mapping_monosomic)

write.table(all_sims, 'tables/power_analysis_complete.tsv', col.names = T, quote = F, row.names = F, sep = '\t')
