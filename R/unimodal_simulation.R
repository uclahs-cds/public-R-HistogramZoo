unimodal.params <- list(
  N = c(25, 500),
  eps = c(0.5, 2),
  noise = c(.05, 0.5)
)

random_unimodal_sim <- function(
    N = c(25, 500),
    unif_length = c(6, 24),
    norm_sd = c(1, 4),
    gamma_shape = c(1, 13),
    eps = c(0.05, 2),
    noise = c(.05, .5),
    max_uniform = NULL,
    remove_low_entropy = NULL,
    truncated_models = FALSE,
    include_data = FALSE,
    include_segment_and_fit = TRUE,
    actual_dist = c('norm', 'unif', 'gamma'),
    metrics = c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq'),
    seed = as.integer(sub('^16', sample(1:9, size = 1), as.integer(Sys.time())))
  ) {

  set.seed(seed);
  sim_dist <- if (length(actual_dist) == 1) actual_dist else sample(actual_dist, size = 1)
  if (is.null(max_uniform)) max_uniform <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(remove_low_entropy)) remove_low_entropy <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(truncated_models)) truncated_models <- sample(c(TRUE, FALSE), size = 1)

  .sample_unif <- function(x, n = 1) {
    runif(n, min = x[1], max = x[2])
  }

  set.seed(seed);
  N_sim <- if (length(N) == 1) N else .sample_unif(N)
  eps_sample <- if (length(eps) == 1) eps else .sample_unif(eps)
  noise_sim <- if (length(noise) == 1) noise else .sample_unif(noise)
  N_noise_sim <- round(N_sim * noise_sim)

  set.seed(seed);
  if (sim_dist == 'norm') {
    param <- if (length(norm_sd) == 1) norm_sd else .sample_unif(norm_sd)
    peak <- rnorm(N_sim, mean = 0, sd = param)
  } else if (sim_dist == 'gamma') {
    param <- if (length(gamma_shape) == 1) gamma_shape else .sample_unif(gamma_shape)
    peak <- rgamma(N_sim, shape = param)
  } else {
    # Unif
    param <-  if (length(unif_length) == 1) unif_length else .sample_unif(unif_length)
    peak <- runif(N_sim, min = 0, max = param)
  }

  data_sd <- sd(peak)
  peak_min <- min(peak)
  peak_max <- max(peak)
  noise_min <- peak_min - data_sd
  noise_max <- peak_max + data_sd

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

  res <- list(
    N = N_sim,
    param = param,
    noise_min = noise_min,
    noise_max = noise_max,
    peak_min = peak_min,
    peak_max = peak_max,
    noise = noise_sim,
    actual_dist = sim_dist,
    eps = eps_sample,
    seed = seed,
    max_uniform = max_uniform,
    remove_low_entropy = remove_low_entropy,
    truncated_models = truncated_models
    )

    print.res <- as.data.frame(res)
    rownames(print.res) <- 'PARAMS: '
    print(print.res)

    all_seg_results <- NULL
    if (include_segment_and_fit) {
      timing <- system.time({
      seg_results_mod <- try({
        segment_and_fit(
          histogram_data,
          max_uniform = max_uniform,
          remove_low_entropy = remove_low_entropy,
          metric = metrics,
          eps = eps_sample,
          truncated_models = truncated_models,
          metric_weights = sqrt(rev(seq_along(metrics))),
          seed = seed
          )
        }, silent = TRUE)
      seg_results <- summarize_results_error(seg_results_mod)
      })

      seg_results <- cbind.data.frame(
        seg_results,
        t(unclass(timing))[, c('user.self', 'sys.self', 'elapsed'), drop = FALSE]
        )

      all_seg_results <- cbind.data.frame(res, seg_results)
      if (!include_data) return(all_seg_results)
      else {
        return(
          list(
            seg_results = all_seg_results,
            hz_model = seg_results_mod,
            peak_data = peak,
            noise_data = noise_data
          )
        )
      }
    }

    if (include_data) {
      # Return the actual data
        res$peak_data <- peak
        res$noise_data <- noise_data
        return(res)
      }
}
