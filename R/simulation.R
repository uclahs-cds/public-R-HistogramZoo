# generate_unimodal_sim_data <- function(

# )

random_unimodal_sim <- function(
    N = c(25, 500),
    unif_length = c(6, 25),
    norm_sd = c(1, 4),
    gamma_shape = c(1, 4),
    eps = c(0.05, 0.5),
    noise = c(.05, .95),
    max_uniform = NULL,
    remove_low_entropy = NULL,
    truncated_models = FALSE,
    metrics = c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq'),
    seed = as.integer(Sys.time())
  ) {
  sim_dist <- sample(c('norm', 'unif', 'gamma'), size = 1)
  if (is.null(max_uniform)) max_uniform <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(remove_low_entropy)) remove_low_entropy <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(truncated_models)) truncated_models <- sample(c(TRUE, FALSE), size = 1)

  .sample_unif <- function(x, n = 1) {
    runif(n, min = x[1], max = x[2])
  }
  set.seed(seed);

  N_sim <- .sample_unif(N)
  eps_sample <- .sample_unif(eps)
  noise_sim <- .sample_unif(noise)
  N_noise_sim <- round(N_sim * noise_sim)

  if (sim_dist == 'norm') {
    param <- .sample_unif(norm_sd)
    peak <- rnorm(N_sim, mean = 0, sd = param)
  } else if (sim_dist == 'gamma') {
    param <- .sample_unif(gamma_shape)
    peak <- rgamma(N_sim, shape = param)
  } else {
    # Unif
    param <- .sample_unif(unif_length)
    peak <- runif(N_sim, min = 0, max = param)
  }

  data_sd <- sd(peak)
  noise_min <- min(peak) - data_sd
  noise_max <- max(peak) + data_sd

  noise_data <- runif(
    N_noise_sim,
    min = noise_min,
    max = noise_max
    )

  peak_noise <- c(
    peak,
    noise_data
    )

  histogram_data <- observations_to_histogram(
    x = peak_noise,
    histogram_bin_width = 1
    )

  timing <- system.time({
    seg_results <- try({
      segment_and_fit(
        histogram_data,
        max_uniform = max_uniform,
        remove_low_entropy = remove_low_entropy,
        metric = metrics,
        eps = eps_sample,
        truncated_models = truncated_models,
        metric_weights = sqrt(rev(seq_along(metrics)))
        )
      }, silent = TRUE)
    seg_results <- summarize_results_error(seg_results)
    })

    seg_results <- cbind.data.frame(
      seg_results,
      t(unclass(timing))[, c('user.self', 'sys.self', 'elapsed'), drop = FALSE]
      )

  res <- list(
    N = N_sim,
    param = param,
    noise_min = noise_min,
    noise_max = noise_max,
    noise = noise_sim,
    actual_dist = sim_dist,
    eps = eps_sample,
    seed = seed,
    max_uniform = max_uniform,
    remove_low_entropy = remove_low_entropy,
    truncated_models = truncated_models
    )

    return(cbind.data.frame(res, seg_results))
}

  # Set optima_flat_endpoints = FALSE
  # Background noise change to M: number of uniform
  # Simulation d: distance between peaks as % of variability

  # Peaks = 1 or 2
  # N = 25... 500
  # Params:
  # Unif: length
  # Norm: sd
  # Gamma: shape (skewness)
  # Distance between peaks d: % overlap of the peaks

  # Alt distribution: Fit semicircle?


  # Simulate epsilon from .05 to .5: function of true distribution and initialization points
  # max_uniform = FALSE
  # remove_low_entropy = TRUE
  # Set default of max_uniform to false?

  # Use the square-root of sqrt(c(5,4,3,2,1)) for weights
  # Determine the best consensus weighting based on the results from the simulation

  # keep consensus method the same
  # Evaluate performance for each metric individually and as consensus
general_sim <- function(
    N,
    param1,
    param2,
    seed = NULL,
    histogram_bin_width = 1,
    rfunc = c('rnorm', 'rgamma', 'runif', 'rlnorm'),
    quiet = TRUE, ...) {
  if (!quiet) {
    cat(sprintf('N = %s\tparam1 = %.3f\tparam2 = %.3f\tseed=%d\n', N, param1, param2, seed))
  }
  rfunc <- get(match.arg(rfunc))

  set.seed(seed)

  histogram_data <- observations_to_histogram(
    rfunc(N, param1, param2),
    histogram_bin_width = histogram_bin_width
    )

  segment_args <- c(
    list(
      histogram_obj = histogram_data,
      seed = seed
      ),
    list(...)
    )

  # Keeping all data and fitting a full model
  timing <- system.time({
    seg_results <- do.call(segment_and_fit, segment_args)
  })

  list(
    timing = timing,
    seg_results = seg_results
  )
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
