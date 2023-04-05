
#' Computes the negative log likelihood of an underlying continuous distribution from binned data.
#'
#' From: https://stats.stackexchange.com/a/68238/97417
#'
#' @param theta parameters for distribution
#' @param cdf the cdf function (pnorm, pgamma, punif)
#' @param counts the histogram bins
#' @param bin_lower the lower bounds for the histogram bins (a, b]
#' @param bin_upper the upper bounds for the histogram bins (a, b]
#' @param ... Additional parameters to pass into `cdf`
#'
#' @return the _negative_ log-likelihood to be _minimized_
#' @export
#'
#' @examples
#' set.seed(100);
#' x <- rnorm(1000, sd = 5)
#' bin.x <- observations_to_histogram(x)
#'
#' bin_log_likelihood(
#'   theta = c(0, 4),
#'   cdf = pnorm,
#'   counts = bin.x$histogram_data,
#'   bin_lower = bin.x$interval_start,
#'   bin_upper = bin.x$interval_end
#'   )
#'
#' dist_optim <- DEoptim::DEoptim(
#'   fn = bin_log_likelihood,
#'   cdf = pnorm,
#'   counts = bin.x$histogram_data,
#'   bin_lower = bin.x$interval_start,
#'   bin_upper = bin.x$interval_end,
#'   lower = c(min(bin.x$interval_start), 0.01),
#'   upper = c(max(bin.x$interval_end), 10),
#'   control = list(
#'     trace = TRUE,
#'     itermax = 500
#'     )
#'   )
#' print(dist_optim$optim$bestmem)
bin_log_likelihood <- function(theta, cdf, counts, bin_lower, bin_upper, ...) {
  common.args <- c(
    as.list(theta),
    list(...)
    )

  res <- unlist(lapply(seq_along(counts), function(i) {
    upper.prob <- do.call(cdf, c(
      list(q = bin_upper[i]),
        common.args
        )
      )
    lower.prob <- do.call(cdf, c(
      list(q = bin_lower[i]),
        common.args
        )
      )
    bin.prob <- upper.prob - lower.prob
    if (!is.na(bin.prob) && bin.prob > 0) counts[i] * log(bin.prob)
    else - Inf
  }))
  - sum(res)
}
