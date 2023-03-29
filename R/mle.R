# From: https://stats.stackexchange.com/a/68238/97417
bin_log_likelihood <- function(theta, cdf, counts, bin_lower, bin_upper) {
  res <- unlist(lapply(seq_along(counts), function(i) {
    bin_prob <- cdf(bin_upper[i], theta[1], theta[2]) - cdf(bin_lower[i], theta[1], theta[2])
    if (!is.na(bin_prob) && bin_prob > 0) counts[i] * log(bin_prob)
    else -Inf
  }))
  - sum(res)
}
