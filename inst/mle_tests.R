library(HistogramZoo)

set.seed(100);
x <- rnorm(1000, sd = 5)
bin.x <- observations_to_histogram(x)
# hist.x <- create_coverageplot(bin.x)
rm(hist.x)

bin_log_likelihood(
  theta = c(0, 4),
  cdf = pnorm,
  counts = bin.x$histogram_data,
  bin_lower = bin.x$interval_start,
  bin_upper = bin.x$interval_end
  )

dist_optim <- DEoptim::DEoptim(
  fn = bin_log_likelihood,
  cdf = pnorm,
  counts = bin.x$histogram_data,
  bin_lower = bin.x$interval_start,
  bin_upper = bin.x$interval_end,
  lower = c(0.1, 0.1),
  upper = c(3, 10),
  control = list(
    trace = TRUE,
    itermax = 500, # Iterations
    VTR = 10^-2 # At 1 %, stop optimizing
    )
  )

x <- rgamma(1000, shape = 3, rate = 0.1)
bin.x <- observations_to_histogram(x)
create_coverageplot(bin.x)

dist_optim_gamma <- DEoptim::DEoptim(
  fn = bin_log_likelihood,
  cdf = pnorm,
  counts = bin.x$histogram_data,
  bin_lower = bin.x$interval_start,
  bin_upper = bin.x$interval_end,
  lower = c(0.1, 0.00001),
  upper = c(5, 3),
  control = list(
    trace = TRUE,
    itermax = 500, # Iterations
    steptol = 50
    )
  )


hist(rtnorm(1e4, mean = 0, sd = 3, a = -1, b = 10))

x <- rgamma(1000, shape = 3, rate = 0.1)
bin.x <- observations_to_histogram(x)

res <- segment_and_fit(
  bin.x,
  distributions = c('gamma', 'norm'),
  metric = 'mle',
  max_uniform = FALSE,
  remove_low_entropy = TRUE,
  truncated_models = TRUE
  )

create_coverageplot(res)
# res$models[[1]]$
