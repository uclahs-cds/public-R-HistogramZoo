
library(HistogramZoo)

# TODO: Add this to testthat

# Parameters
metrics <- c('jaccard', 'intersection', 'ks', 'mse', 'chisq')

unimodal_params <- list(
  N = c(25, 500),
  unif_length = c(6, 24),
  norm_sd = c(1, 4),
  gamma_shape = c(1, 13),
  eps = c(0.5, 2),
  noise = c(.05, .5),
  max_uniform = NULL,
  remove_low_entropy = NULL,
  truncated_models = FALSE,
  metrics = metrics,
  run_segment_and_fit = TRUE,
  include_simulated_data = FALSE,
  include_fit_data = TRUE
)

multimodal_params <- c(
  unimodal_params,
  list(
    peaks = 2:4,
    peak_shift = c(2, 8) # Peak shift in standard deviations from previous peak
  )
)

unimodal_seed <- 343726688
unimodal_params <- c(unimodal_params, list(seed = unimodal_seed))
multimodal_seed <- 645469583
multimodal_params <- c(multimodal_params, list(seed = multimodal_seed))

uni <- do.call(HistogramZoo:::random_unimodal_sim, unimodal_params)
multi <- do.call(HistogramZoo:::random_multi_peak_sim, multimodal_params)