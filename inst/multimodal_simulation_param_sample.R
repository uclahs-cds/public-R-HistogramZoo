library(HistogramZoo)

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results', 'multi_modal_sim_v2');
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

params <- list(
    N = c(25, 500),
    unif_length = c(6, 25),
    norm_sd = c(1, 4),
    gamma_shape = c(1, 4),
    eps = c(0.5, 2),
    noise = c(.05, 0.5),
    max_uniform = NULL,
    remove_low_entropy = NULL,
    truncated_models = FALSE,
    peaks = 2:4,
    peak_shift = c(1, 5), # Peak shift in standard deviations from previous peak
    metrics = c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq')
    )
cat('Params: \n')
print(params)

sim.data <- replicate(
  N,
  expr = {
    cat('V2 New multi-modal sim: ', as.character(Sys.time()), '\n');
    tryCatch(
      expr = {
        do.call( HistogramZoo:::random_multi_peak_sim, params);
      },
      error = function(x) {
        cat(x);
        }
      )
    },
  simplify = FALSE
  )

sim.data <- sim.data[sapply(sim.data, Negate(is.null))];

suffix <- '_multi-peak-sim_'
if (mle) suffix <- paste0(suffix, '_MLE_')
file.prefix <- gsub('[ ]', '_', as.integer(Sys.time()))

results.names <- names(sim.data[[1]])
for (n in results.names) {
  filename <- paste0(file.prefix, suffix, n, '.tsv')
  result.df <- do.call(
    plyr::rbind.fill,
    lapply(sim.data, function(x) {
      x[[n]]
    })
  )

  result.df$mle <- mle

  full.filename <- file.path(results.folder, filename)
  cat('Saving output to: ', full.filename, '\n')
  write.table(
    x  = result.df,
    file = full.filename,
    sep = '\t',
    row.names = FALSE
    )
  }
