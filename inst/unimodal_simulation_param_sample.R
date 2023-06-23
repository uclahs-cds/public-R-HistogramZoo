library(HistogramZoo)

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results', 'unimodal_sim_noise_mle');
set.seed(314)

opts <- commandArgs(trailingOnly = TRUE)
N <- 10
if (length(opts) >= 1) {
  N <- opts[[1]]
}

sim.data <- replicate(
  N,
  HistogramZoo:::random_unimodal_sim(
    metric = 'mle'
  ),
  simplify = FALSE
  )

sim.df <- do.call(
    plyr::rbind.fill,
    lapply(seq_along(sim.data), function(j) {
        HistogramZoo:::process_sim(sim.data[[j]])
        })
    )

filename <- gsub('[ ]', '_', paste0(Sys.time(), '_', sample(1:1e4, 1), '_Unimodal_Sim.tsv'))
full.filename <- file.path(results.folder, filename)
cat('Saving output to: ', full.filename, '\n')
# saveRDS(sim.df, file = full.filename)
write.table(sim.df, file = full.filename, sep = '\t', row.names = FALSE)
