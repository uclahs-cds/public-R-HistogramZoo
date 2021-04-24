fit.continuous.distributions = function(
  x,
  sd.scale,
  seg.start,
  seg.end,
  fit.normal.mixture = T,
  max.iterations = 500
){

  # Initializing
  mod = list()

  # Uniform Distribution
  tryCatch(
    expr = {
      mod$unif <- fitdistrplus::fitdist(
        data = x,
        distr = "unif")
      mod$unif$sd_scale <- sd.scale
    },
    error = function(e) {
      warning(sprintf("Error in fitdist unif for segment [%d, %d]", seg.start, seg.end))
      print(e)
    }
  )

  # Truncated Normal Distribution
  tryCatch(
    expr = {
      mod$tnorm <- fitdistrplus::fitdist(
        data = x,
        distr = "tnorm",
        fix.arg = list(a = 0, b = max(x) + 1e-10),
        start = list(mean = mean(x), sd = sd(x)),
        optim.method="L-BFGS-B")
      mod$tnorm$sd_scale <- sd.scale
    },
    error = function(e) {
      warning(sprintf("Error in fitdist tnorm for segment [%d, %d]", seg.start, seg.end))
      print(e)
    }
  )

  # Truncated gamma distribution
  tryCatch(
    expr = {
      mod$tgamma <- fitdistrplus::fitdist(
        data = x,
        distr = "tgamma",
        fix.arg = list(a = 0, b = max(x)),
        start = list(shape = 2, rate = 1),
        # Set the lower bound for shape and rate params
        lower = c(1, 0.5))
      mod$tgamma$sd_scale <- sd.scale
    },
    error = function(e) {
      warning(sprintf("Error in fitdist tgamma for segment [%d, %d]", seg.start, seg.end))
      print(e)
    }
  )

  # Flipped Truncated gamma distribution
  tryCatch(
    expr = {
      mod$tgamma_flip <- fitdistrplus::fitdist(
        data = x,
        distr = "tgamma_flip",
        fix.arg = list(b = max(x) + 1e-10),
        start = list(shape = 2, rate = 1),
        # Set the lower bound for shape and rate params
        lower = c(1, 0.5))
      mod$tgamma_flip$sd_scale <- sd.scale
    },
    error = function(e) {
      warning(sprintf("Error in fitdist tgamma_flip for segment [%d, %d]", seg.start, seg.end))
      print(e)
    }
  )

  # Mixture of Normals
  if(fit.normal.mixture) {
    # Fit a Normal Mixture model
    maxiter <- max.iterations
    # Fit mixture model, silencing output
    invisible(capture.output({
      mixfit <- mixtools::normalmixEM(x, verb = FALSE, maxit = maxiter, epsilon = 1e-04, k = 2)
    }))
    if((length(mixfit$all.loglik) - 1) >= maxiter) {
      # EM did not converge. Don't use results.
      mixfit <- NULL
    } else {
      mod$norm_mixture <- mixfit
      mod$norm_mixture$sd_scale <- sd.scale
    }
  }

  return(mod)

}
