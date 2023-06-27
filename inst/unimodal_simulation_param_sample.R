library(HistogramZoo)

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results', 'unimodal_sim_noise_v2');
set.seed(314)

opts <- commandArgs(trailingOnly = TRUE)
N <- 10
if (length(opts) >= 1) {
  N <- opts[[1]]
}

metrics <- c('jaccard', 'intersection', 'ks', 'mse', 'chisq')

mle <- FALSE
if (length(opts) >= 2) {
    mle <- tolower(opts[[2]]) == 'mle'
    cat('Using MLE\n');
    metrics <- 'mle'
    }

sim.data <- replicate(
  N,
  {
    cat('New unimodal sim with : ', as.character(Sys.time()), '\n');
    HistogramZoo:::random_unimodal_sim(
        metrics = metrics
        )
  },
  simplify = FALSE
  )

sim.df <- do.call(
    plyr::rbind.fill,
    sim.data
    )


suffix <- '_Unimodal_Sim'
if (mle) suffix <- paste0(suffix, '_MLE')

filename <- gsub('[ ]', '_', paste0(Sys.time(), '_', sample(1:1e4, 1), suffix, '.tsv'))
full.filename <- file.path(results.folder, filename)
cat('Saving output to: ', full.filename, '\n')
write.table(sim.df, file = full.filename, sep = '\t', row.names = FALSE)
