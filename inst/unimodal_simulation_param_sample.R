library(HistogramZoo)

output_path <- '/hot/users/stefaneng/public-R-HistogramZoo/results'
set.seed(314)

opts <- commandArgs(trailingOnly = TRUE)
N <- 10
if (length(opts) > 1) {
  N <- opts[[1]]
}


results <- replicate(
  N,
  random_unimodal_sim(),
  simplify = FALSE
  )

filename <- paste0(Sys.time(), '_', sample(1:1e4, 1), '_Unimodal_Sim.rds')
saveRDS(resutls, file = file.path(output_path, filename))
