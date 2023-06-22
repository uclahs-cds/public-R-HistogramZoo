library(HistogramZoo)

output_path <- '/hot/users/stefaneng/public-R-HistogramZoo/results/unimodal_sim_noise'
set.seed(314)

opts <- commandArgs(trailingOnly = TRUE)
N <- 10
if (length(opts) >= 1) {
  N <- opts[[1]]
}

sim.data <- replicate(
  N,
  HistogramZoo:::random_unimodal_sim(),
  simplify = FALSE
  )

process_sim <- function(x) {
    res <- x[c('N', 'dist', 'param', 'noise', 'eps')]
    names(res) <- c('N', 'actual_dist', 'param', 'noise', 'eps')
    res$timing <- x$timing[['user.self']]
    if (! is(x$seg_results,  'try-error')) {
        models <- x$seg_results$models
        metrics <- x$seg_results$metric
        res$num_segments <- length(models)
        all_metric_results <- do.call(plyr::rbind.fill, lapply(metrics, summarize_results, result = x$seg_results))

        cbind.data.frame(all_metric_results, res)
    } else {
        cbind.data.frame( data.frame(error = x[[1]]), res)
    }
}

sim.df <- do.call(
    plyr::rbind.fill,
    lapply(seq_along(sim.data), function(j) {
        process_sim(sim.data[[j]])
        })
    )

filename <- gsub('[ ]', '_', paste0(Sys.time(), '_', sample(1:1e4, 1), '_Unimodal_Sim.tsv'))
full.filename <- file.path(output_path, filename)
cat('Saving output to: ', full.filename, '\n')
# saveRDS(sim.df, file = full.filename)
write.table(sim.df, file = full.filename, sep = '\t', row.names = FALSE)
