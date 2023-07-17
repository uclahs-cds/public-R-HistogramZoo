library(HistogramZoo)

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results', 'unimodal_sim_noise_v4');
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

params <- list(
    N = c(25, 500),
    unif_length = c(6, 25),
    norm_sd = c(1, 4),
    gamma_shape = c(1, 4),
    eps = c(0.5, 2),
    noise = c(.05, .5),
    max_uniform = NULL,
    remove_low_entropy = NULL,
    truncated_models = FALSE,
    metrics = metrics
    )

cat('Params: \t')
print(params)

sim.data <- replicate(
  N,
  {
    cat('New unimodal sim: ', as.character(Sys.time()), '\n');
    do.call(
      HistogramZoo:::random_unimodal_sim,
      params
      )
  },
  simplify = FALSE
  )

sim.df <- do.call(
    plyr::rbind.fill,
    sim.data
    )


suffix <- 'Unimodal_Sim'

if (mle) suffix <- paste0(suffix, '_MLE')
file.prefix <- gsub('[ ]', '_', as.integer(Sys.time()))

filename <- gsub('[ ]', '_', paste0(file.prefix, '_', suffix, '.tsv'))
full.filename <- file.path(results.folder, filename)
cat('Saving output to: ', full.filename, '\n')
write.table(sim.df, file = full.filename, sep = '\t', row.names = FALSE)
