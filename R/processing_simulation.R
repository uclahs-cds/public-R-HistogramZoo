# This is just temporary to add peak_min/peak_max to missing sims
peak_min_recompute <- function(
    N = c(25, 500),
    unif_length = c(6, 25),
    norm_sd = c(1, 4),
    gamma_shape = c(1, 4),
    eps = c(0.05, .5),
    noise = c(.05, .5),
    max_uniform = c(TRUE, FALSE),
    remove_low_entropy = c(TRUE, FALSE),
    truncated_models = c(TRUE, FALSE),
    actual_dist = c('norm', 'unif', 'gamma'),
    metrics = c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq'),
    seed = as.integer(Sys.time()),
    ...
  ) {

  .sample_unif <- function(x, n = 1) {
    runif(n, min = x[1], max = x[2])
  }
  set.seed(seed);

  N_sim <- .sample_unif(N)
  eps_sample <- .sample_unif(eps)
  noise_sim <- .sample_unif(noise)
  N_noise_sim <- round(N_sim * noise_sim)

  if (actual_dist == 'norm') {
    param <- .sample_unif(norm_sd)
    peak <- rnorm(N_sim, mean = 0, sd = param)
  } else if (actual_dist == 'gamma') {
    param <- .sample_unif(gamma_shape)
    peak <- rgamma(N_sim, shape = param)
  } else {
    # Unif
    param <- .sample_unif(unif_length)
    peak <- runif(N_sim, min = 0, max = param)
  }

  data_sd <- sd(peak)
  peak_min <- min(peak)
  peak_max <- max(peak)
  noise_min <- peak_min - data_sd
  noise_max <- peak_max + data_sd
  list(seed = seed, peak_min = peak_min, peak_max = peak_max, noise_min = noise_min, noise_max = noise_max)
}

process_sim <- function(x) {
    res <- x[c('N', 'dist', 'param', 'noise', 'eps')]
    names(res) <- c('N', 'actual_dist', 'param', 'noise', 'eps')
    res$timing <- x$timing[['user.self']]

    if (class(x$seg_results) == 'HistogramFit') {
      seg.results <- summarize_results_error(x$seg_results)
    } else {
      seg.results <- do.call(plyr::rbind.fill, x$seg_results)
    }

    return(cbind.data.frame(res, seg.results))
}

summarize_results_error <- function(x) {
  if (! is(x,  'try-error')) {
        models <- x$models
        metrics <- x$metric
        all_metric_results <- do.call(plyr::rbind.fill, lapply(metrics, summarize_results, result = x))

        rtn <- all_metric_results
        rtn$num_segments <- length(models)
        rtn
    } else {
        rtn <- data.frame(error = x[[1]])
    }
  return(rtn);
}

overlap_size <- function(a1, a2, b1, b2) {
  max(0, min(a2, b2) - max(a1, b1))
}

union_size <- function(a1, a2, b1, b2) {
  if (overlap_size(a1, a2, b1, b2) > 0) max(a2, b2) - min(a1, b1)
  else (a2 - a1) + (b2 - b1)
}
