general_sim <- function(
    N,
    param1,
    param2,
    seed = 131313,
    rfunc = c('rnorm', 'rgamma', 'runif', 'rlnorm'),
    quiet = TRUE, ...) {
  if (!quiet) {
    cat(sprintf('N = %s\tparam1 = %.3f\tparam2 = %.3f\tseed=%d\n', N, param1, param2, seed))
  }
  rfunc <- get(match.arg(rfunc))

  set.seed(seed)

  histogram_data <- observations_to_histogram(
    rfunc(N, param1, param2)
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
