library(HistogramZoo)

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results', 'multi_modal_sim');
set.seed(314)

opts <- commandArgs(trailingOnly = TRUE)
N <- 100
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
  expr = {
    cat('New multi-modal sim: ', as.character(Sys.time()), '\n');
    HistogramZoo:::random_multi_peak_sim(
      metrics = metrics,
      peaks = 2:4,
      return_fit = FALSE
      )
    },
  simplify = FALSE
  )

suffix <- '_multi-peak-sim_'
if (mle) suffix <- paste0(suffix, '_MLE')
file.prefix <- gsub('[ ]', '_', as.integer(Sys.time()))

results.names <- names(sim.data[[1]])
for (n in results.names) {
  filename <- paste0(file.prefix, suffix, n, '.tsv')
  result.df <- do.call(
    'rbind.data.frame',
    lapply(sim.data, function(x) {
      x[[n]]
    })
  )

  full.filename <- file.path(results.folder, filename)
  cat('Saving output to: ', full.filename, '\n')
  write.table(
    x  = result.df,
    file = full.filename,
    sep = '\t',
    row.names = FALSE
    )
  }
