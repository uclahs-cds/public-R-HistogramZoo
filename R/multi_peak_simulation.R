random_multi_peak_sim <- function(
    N = c(25, 500),
    unif_length = c(6, 25),
    norm_sd = c(1, 4),
    gamma_shape = c(1, 4),
    eps = c(0.05, 0.5),
    noise = c(.05, .5),
    max_uniform = NULL,
    remove_low_entropy = NULL,
    truncated_models = FALSE,
    peaks = 2:3,
    peak_shift = c(1, 5), # Peak shift in standard deviations from previous peak
    metrics = c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq'),
    seed = as.integer(Sys.time())
  ) {
  if (is.null(max_uniform)) max_uniform <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(remove_low_entropy)) remove_low_entropy <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(truncated_models)) truncated_models <- sample(c(TRUE, FALSE), size = 1)

  .sample_unif <- function(x, n = 1) {
    runif(n, min = x[1], max = x[2])
  }
  set.seed(seed);

  if (length(peaks) > 1) {
    peaks_sim <- sample(x = peaks, size = 1)
    } else {
    peaks_sim <- peaks
    }

  sim_dist <- sample(c('norm', 'unif', 'gamma'), size = peaks_sim, replace = TRUE)
  N_sim <- round(.sample_unif(N, n = peaks_sim))
  eps_sample <- .sample_unif(eps, n = 1)
  noise_sim <- .sample_unif(noise)
  N_noise_sim <- round(sum(N_sim) * noise_sim)

  # Shift...

  multi_histogram_data <- lapply(seq_len(peaks_sim), function(i) {
    # Mean is 0 for each of the distributions
    if (sim_dist[[i]] == 'norm') {
      param <- .sample_unif(norm_sd)
      peak <- rnorm(N_sim[[i]], mean = 0, sd = param)
      sd_dist <- param
    } else if (sim_dist[[i]] == 'gamma') {
      param <- .sample_unif(gamma_shape)
      peak <- rgamma(N_sim[[i]], shape = param)
      # Subtract of param = shape = mean for gamma distribution
      peak <- peak - param
      sd_dist <- sqrt(param) # When rate = 1
    } else {
      # Unif
      param <- .sample_unif(unif_length)
      peak <- runif(N_sim[[i]], min = -param/2, max = param/2)
      sd_dist <- param * (1 / sqrt(12))
    }

    # To get information about peak
    # min, max, etc.
    peak_hist <- observations_to_histogram(
      peak,
      histogram_bin_width = 1
      )

    list(
      param = param,
      peak = peak,
      sd_sample = sd(peak),
      sd_dist = sd_dist
      )
    })

  peak_shift_sample <- 0
  # Shift each eap
  # In standard deviations
  if (peaks_sim > 1) {
    peak_shift_sample <- c(0, .sample_unif(peak_shift, n = peaks_sim - 1))
    for (p in 2:peaks) {
      # Shift peaks
      multi_histogram_data[[p]]$shift <- peak_shift_sample[p]
      multi_histogram_data[[p]]$peak <- multi_histogram_data[[p]]$peak +
        peak_shift_sample[p] * multi_histogram_data[[p]]$sd_dist
    }
  }

  # Make histograms from peaks

  # Raw un-binned peak sd
  merged_raw_data <- unlist(sapply(multi_histogram_data, function(x) x$peak))
  sd_raw_data <- sd(merged_raw_data)

  peak_histograms <- lapply(multi_histogram_data, function(x) {
    observations_to_histogram(x$peak)
    })

  # TODO: Can we assign % of each peak into each bin to assess how well segmentation worked?
  # lapply(peak_histograms, function(x) x$interval_start)

  merged_peak_histograms <- Reduce(`+`, peak_histograms)

  noise_min <- min(merged_raw_data) - sd_raw_data
  noise_max <- max(merged_raw_data) + sd_raw_data

  noise_data <- runif(
      N_noise_sim,
      min = noise_min,
      max = noise_max
      )

  noise_histogram <- observations_to_histogram(noise_data)

  # merged_peaks <- unlist(sapply(multi_histogram_data, function(x) x$peak))
  merged_peak_histograms_noise <- merged_peak_histograms +
    observations_to_histogram(
        x = noise_data,
        histogram_bin_width = 1
        )

  timing <- system.time({
    seg_results <- try({
      segment_and_fit(
        merged_peak_histograms_noise,
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

  peak_min <- unlist(lapply(peak_histograms, function(x) head(x$interval_start, n = 1)))
  peak_max <- unlist(lapply(peak_histograms, function(x) tail(x$interval_end, n = 1)))

  res <- list(
    N = N_sim,
    param = unlist(lapply(multi_histogram_data, function(x) x$param)),
    peak_min = peak_min,
    peak_max = peak_max,
    peak_shift = peak_shift_sample,
    noise_min = noise_min,
    noise_max = noise_max,
    noise = noise_sim,
    actual_dist = sim_dist,
    eps = eps_sample,
    seed = seed,
    max_uniform = max_uniform,
    remove_low_entropy = remove_low_entropy,
    truncated_models = truncated_models,
    peaks = peaks_sim
    )

    # long format, one peak per line
    res <- data.frame(res)

    # Joined by seed
    seg_results$seed <- seed

    return(
      list(
        actual_peaks = res,
        seg_results = seg_results
        )
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
