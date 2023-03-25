library(HistogramZoo)

output_path <- '/hot/users/stefaneng/public-R-HistogramZoo/results'
set.seed(314)

opts <- commandArgs(trailingOnly = TRUE)
N <- 10
if (length(opts) >= 1) {
  N <- opts[[1]]
}


results <- replicate(
  N,
  HistogramZoo:::random_unimodal_sim(),
  simplify = FALSE
  )

filename <- gsub('[ ]', '_', paste0(Sys.time(), '_', sample(1:1e4, 1), '_Unimodal_Sim.rds'))
full_filename <- file.path(output_path, filename)
cat('Saving output to: ', full_filename, '\n')
saveRDS(results, file = full_filename)
