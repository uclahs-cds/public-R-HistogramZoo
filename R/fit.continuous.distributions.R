fit.continuous.distributions = function(
  x,
  seg.start,
  seg.end,
  fit.mixtures = c("unif", "tnorm", "tgamma", "tgamma_flip", "mixEM"),
  max.iterations = 500
){

  # Calculating initializing values for parameter estimation
  x.mean = mean(x)
  x.var = var(x)
  x.sd = sd(x)
  x.shape.initial = x.mean^2/x.var
  x.rate.initial = x.mean/x.var

  # Initializing list of models
  mod = list()

  # Uniform Distribution
  if("unif" %in% fit.mixtures){
    tryCatch(
      expr = {
        mod$unif <- fitdistrplus::fitdist(
          data = x,
          distr = "unif")
      },
      error = function(e) {
        warning(sprintf("Error in fitdist unif for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )
  }

  # Truncated Normal Distribution
  if("tnorm" %in% fit.mixtures){
    tryCatch(
      expr = {
        mod$tnorm <- fitdistrplus::fitdist(
          data = x,
          distr = "tnorm",
          fix.arg = list(a = 0, b = max(x) + 1e-10),
          start = list(mean = x.mean, sd = x.sd),
          optim.method="L-BFGS-B")
      },
      error = function(e) {
        warning(sprintf("Error in fitdist tnorm for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )
  }

  # Truncated gamma distribution
  if("tgamma" %in% fit.mixtures){
    tryCatch(
      expr = {
        mod$tgamma <- fitdistrplus::fitdist(
          data = x,
          distr = "tgamma",
          fix.arg = list(a = 0, b = max(x)),
          start = list(shape = x.shape.initial, rate = x.rate.initial),
          # Set the lower bound for shape and rate params
          lower = c(1, 0.00001))
      },
      error = function(e) {
        warning(sprintf("Error in fitdist tgamma for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )
  }

  # Flipped Truncated gamma distribution
  if("tgamma_flip" %in% fit.mixtures){
    tryCatch(
      expr = {
        mod$tgamma_flip <- fitdistrplus::fitdist(
          data = x,
          distr = "tgamma_flip",
          fix.arg = list(b = max(x) + 1e-10),
          start = list(shape = x.shape.initial, rate = x.rate.initial),
          # Set the lower bound for shape and rate params
          lower = c(1, 0.00001))
      },
      error = function(e) {
        warning(sprintf("Error in fitdist tgamma_flip for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )
  }

  # Mixture of Normals
  if("mixEM" %in% fit.mixtures) {
    # Fit a Normal Mixture model
    maxiter <- max.iterations
    # Fit mixture model, silencing output
    out = tryCatch(
      {
      invisible(capture.output({
        mixfit <- mixtools::normalmixEM(x, verb = FALSE, maxit = maxiter, epsilon = 1e-04, k = 2)
      }))
      mixfit
      },
      error = function(e) {
        print(e)
        return(NULL)
      }
    )
    mod$norm_mixture <- out
  }

  return(mod)

}

# Fit the model parameters by optimizing a histogram metric
fit.distributions.optim <- function(x, metric = c("jaccard", "intersection", "ks"), truncated = FALSE) {
  metric = match.arg(metric)
  # Get one of the metrics from histogram.distances
  metric.func = get(paste('histogram', metric, sep = "."))
  bin = x[, 1]
  freq = x[, 2]
  N = sum(freq)

  hist.optim <- function(params, dist = c("norm", "gamma")) {
  # Compute the expected counts for the given parameters
    dist = match.arg(dist)
    args = c(list(x = bin), params)
    if(truncated) {
      args$a =  min(bin) - 1e-10
      args$b = max(bin) + 1e-10
    }
    trunc.letter = if(truncated) "t" else ""
    dens = do.call(paste0("d", trunc.letter, dist), args) * N
    dens[is.na(dens)] = 0
    res = metric.func(freq, dens)
    if(is.na(res) || res == -Inf) browser()
    res
  }

  L = sum(bin)
  hist.mean = sum(freq * bin) / L
  hist.var = sum(freq * (bin - hist.mean)^2) / L
  norm.res <- optim(par = c(hist.mean, sqrt(hist.var)),
                    fn = hist.optim,
                    method = "L-BFGS-B",
                    dist = "norm",
                    # fnscale = -1 does maximization and 1 for minimization
                    control = list(fnscale = 1),
                    lower = c(-Inf, 0.001),
                    upper = c(max(bin), (max(bin) - min(bin)) * 0.5))

  shape.init = hist.mean^2 / hist.var
  rate.init = hist.mean / hist.var
  gamma.res <- optim(par = c(shape.init, rate.init),
                    fn = hist.optim,
                    method = "L-BFGS-B",
                    dist = "gamma",
                    control = list(fnscale = 1),
                    lower = c(0.001, 0.01))
  list(tnorm = norm.res,
       tgamma = gamma.res)
}
