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
