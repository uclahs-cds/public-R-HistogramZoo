library(HistogramZoo);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results', 'unimodal_sim_noise');
sim.rds.paths <- list.files(results.folder, pattern = 'Unimodal_Sim.rds$', full.names = TRUE);

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


all_results <- do.call(
    plyr::rbind.fill,
    lapply(seq_along(sim.rds.paths), function(i) {
        cat('[', i, '] ', 'Processing: ', sim.rds.paths[[i]], '\n')
        sim.data <- readRDS(sim.rds.paths[[i]])
        sim.df <- do.call(
            plyr::rbind.fill,
            lapply(seq_along(sim.data), function(j) {
                res <- process_sim(sim.data[[j]])
                res$id <- (i - 1) * 10 + j
                res
                })
            )
            write.table(
                x = sim.df,
                file = paste0(sim.rds.paths[[i]], '_processed_results.tsv'),
                sep = '\t',
                row.names = FALSE
                )
            sim.df
            })
    )

# TODO: Split this into batches and combine?
write.table(
    x = all_results,
    file = file.path(results.folder, 'processed_results.tsv'),
    sep = '\t',
    row.names = FALSE
    )
