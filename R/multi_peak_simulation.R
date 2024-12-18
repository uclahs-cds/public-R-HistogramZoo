#' Simulate a multimodal histogram composed as the sum of data from sampling several distributions
#' WARNING: Comprehensive error checking required before export
#' TODO: potentially allow for user specification of sampled distributions (i.e. actual dist)
#'
#' @param N number of points to sample from the distribution
#' @param unif_length range of uniform distribution length (vector of length 2)
#' @param norm_sd range of normal standard deviations (vector of length 2)
#' @param gamma_shape range of gamma shape (vector of length 2)
#' @param eps range of hyperparameter epsilon (vector of length 2)
#' @param noise range of noise as a proportion of total counts (vector of length 2)
#' @param max_uniform NULL forces an equal probability of sampling TRUE or FALSE, otherwise deterministic
#' @param remove_low_entropy NULL forces an equal probability of sampling TRUE or FALSE, otherwise deterministic
#' @param truncated_models NULL forces an equal probability of sampling TRUE or FALSE, otherwise deterministic
#' @param peaks range of number of peaks to sample from (vector of length 2)
#' @param peak_shift range of number of standard deviations of the first peak to sample from (vector of length 2) - this is used to shift subsequent peaks
#' @param metrics metrics for segment and fit
#' @param seed numeric seed
#' @param run_segment_and_fit logical; whether to run segment and fit or not
#' @param include_simulated_data logical; whether to return simulated data
#' @param include_fit_data logical; whether to return the HistogramFit object
#'
#' @return A list of results
#' \describe{
#'     \item{actual_peaks}{a data frame of parameters used for simulation, one row per simulated peak}
#'     \item{seg_results}{if `run_segment_and_fit`, a data frame for peaks identified by segment_and_fit}
#'     \item{overlap}{if `run_segment_and_fit`, a data frame for overlap between simulated and fitted peaks}
#'     \item{peak_data}{if `include_simulated_data`, simulated peak counts will be returned}
#'     \item{noise_data}{if `include_simulated_data`, simulated noise counts will be returned}
#'     \item{hz_model}{if `include_fit_data`, the HistogramFit object will be returned}
#'}
#'
#' @examples \dontrun{
#' simulated_histogram <- random_multi_peak_sim()
#'}
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
    peak_shift = c(2, 8), # Peak shift in standard deviations from first peak
    metrics = c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq'),
    seed = as.integer(sub('^16', sample(1:9, size = 1), as.integer(Sys.time()))),
    run_segment_and_fit = TRUE,
    include_simulated_data = FALSE,
    include_fit_data = FALSE
  ) {

  set.seed(seed);
  if (is.null(max_uniform)) max_uniform <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(remove_low_entropy)) remove_low_entropy <- sample(c(TRUE, FALSE), size = 1)
  if (is.null(truncated_models)) truncated_models <- sample(c(TRUE, FALSE), size = 1)

  .sample_unif <- function(x, n = 1) {
    runif(n, min = x[1], max = x[2])
  }

  if (length(peaks) > 1) {
    peaks_sim <- sample(x = peaks, size = 1)
    } else {
    peaks_sim <- peaks
    }

  set.seed(seed);
  sim_dist <- sample(c('norm', 'unif', 'gamma'), size = peaks_sim, replace = TRUE)
  N_sim <- round(.sample_unif(N, n = peaks_sim))
  eps_sample <- .sample_unif(eps, n = 1)
  noise_sim <- .sample_unif(noise)
  N_noise_sim <- round(sum(N_sim) * noise_sim)

  set.seed(seed);
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
  # Shift each peak (a sampled number of standard deviations)
  if (peaks_sim > 1) {
    peak_shift_sample <- cumsum(c(0, .sample_unif(peak_shift, n = peaks_sim - 1)))
    peaK_shift_actual <- multi_histogram_data[[1]]$sd_dist * peak_shift_sample
    for (p in 2:peaks_sim) {
      # Shift peaks
      multi_histogram_data[[p]]$shift <- peak_shift_sample[p]
      multi_histogram_data[[p]]$true_shift <- peaK_shift_actual[p]
      multi_histogram_data[[p]]$peak <- multi_histogram_data[[p]]$peak +
        multi_histogram_data[[p]]$true_shift
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

  peak_min <- unlist(lapply(peak_histograms, function(x) head(x$interval_start, n = 1)))
  peak_max <- unlist(lapply(peak_histograms, function(x) tail(x$interval_end, n = 1)))

  optima <- find_local_optima(merged_peak_histograms_noise)

    res <- list(
        N = N_sim,
        local_optima = sum(sapply(optima, length)),
        param = unlist(lapply(multi_histogram_data, function(x) x$param)),
        peak_num = seq_along(peak_min),
        peak_min = peak_min,
        peak_max = peak_max,
        peak_shift = peak_shift_sample,
        peak_shift_actual = peaK_shift_actual,
        noise_min = noise_min,
        noise_max = noise_max,
        noise = noise_sim,
        noise_N = N_noise_sim,
        actual_dist = sim_dist,
        eps = eps_sample,
        seed = seed,
        max_uniform = max_uniform,
        remove_low_entropy = remove_low_entropy,
        truncated_models = truncated_models,
        peaks = peaks_sim
        )

    # long format, one peak per line
    actual_peaks <- as.data.frame(res)

    print('PARAMS: ')
    print(actual_peaks)

  if (run_segment_and_fit) {
    timing <- system.time({
      seg_results_mod <- try({
        segment_and_fit(
          merged_peak_histograms_noise,
          max_uniform = max_uniform,
          remove_low_entropy = remove_low_entropy,
          metric = metrics,
          eps = eps_sample,
          truncated_models = truncated_models,
          metric_weights = sqrt(rev(seq_along(metrics))),
          seed = seed
        )
      }, silent = TRUE)
    peak_start <- seg_results_mod$interval_start[seg_results_mod$p[, 'start']]
    peak_end <- seg_results_mod$interval_end[seg_results_mod$p[, 'end']]
    seg_results <- summarize_results_error(seg_results_mod)
    })

    seg_results <- cbind.data.frame(
      seg_results,
      t(unclass(timing))[, c('user.self', 'sys.self', 'elapsed'), drop = FALSE]
    )

    # Joined by seed
    seg_results$seed <- seed

    # Compute overlap data for peaks
    overlap_data <- merge(
      actual_peaks[, c('peak_num', 'peak_min', 'peak_max')],
      data.frame(
        fit_peak_start = peak_start,
        fit_peak_end = peak_end,
        fit_peak_num = seq_along(peak_start)
        ),
      all = TRUE,
      by = character() # Cross join
      )

    overlap_data$overlap <- apply(
      X = overlap_data[, c('peak_min', 'peak_max', 'fit_peak_start', 'fit_peak_end')],
      MARGIN = 1,
      FUN = function(x) {
        do.call(overlap_size, as.list(unname(x)))
      })

    overlap_data <- data.frame(
      overlap_data,
      actual_dist = sim_dist,
      eps = eps_sample,
      param = actual_peaks$param,
      N = N_sim,
      noise = noise_sim,
      noise_N = N_noise_sim,
      max_uniform = max_uniform,
      remove_low_entropy = remove_low_entropy,
      truncated_models = truncated_models,
      seed = seed
      )
  }

  rtn <- list(
    "actual_peaks" = actual_peaks
  )

  if(run_segment_and_fit){
    rtn <- c(
      rtn,
      list(
        "seg_results" = seg_results,
        "overlap" = overlap_data
      )
    )
  }

  if(include_simulated_data){
    rtn <- c(
      rtn,
      list(
        "peak_data" = merged_raw_data,
        "noise_data" = noise_data
      )
    )
  }

  if(include_fit_data && run_segment_and_fit){
    rtn <- c(
      rtn,
      list(
        "hz_model" = seg_results_mod
      )
    )
  }

  return(rtn)
}
