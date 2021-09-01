#' Fit the model parameters by optimizing a histogram metric
fit.distributions.optim_old = function(x, metric = c("jaccard", "intersection", "ks", "mse", "chisq"), truncated = FALSE, distr = c("norm", "gamma", "unif")) {
  metric = match.arg(metric)
  distr = match.arg(distr, several.ok = TRUE)
  # Get one of the metrics from histogram.distances
  metric.func = get(paste('histogram', metric, sep = "."))
  bin = x[, 1]
  freq = x[, 2]
  N = sum(freq)
  L = sum(bin)
  hist.mean = sum(freq * bin) / N
  hist.var = sum(freq * (bin - hist.mean)^2) / N

  .hist.optim = function(params, .dist = c("norm", "gamma")) {
    # Compute the expected counts for the given parameters
    args = c(list(x = bin), params)
    if(truncated) {
      args$a =  min(bin) - 1e-10
      args$b = max(bin) + 1e-10
    }
    trunc.letter = if(truncated) "t" else ""
    dens = tryCatch({
      do.call(paste0("d", trunc.letter, .dist), args) * N
    }, error = function(err) {
      rep(0, length(bin))
    })
    dens[is.na(dens)] = 0
    res = metric.func(freq, dens)
    if(is.na(res) || res == -Inf || res == Inf) browser()
    res
  }

  rtn = list()
  if("unif" %in% distr) {
    # Add uniform distribution
    unif.dens = 1 / (max(bin) - min(bin))
    rtn$unif = list(
      dist = "unif",
      value = metric.func(freq, rep(unif.dens * N, length(L))),
      dens = function(x = NULL, scale = TRUE) {
        if(missing(x)) {
          x = bin
        }
        res = ifelse(x >= min(bin) & x <= max(bin), unif.dens, 0)
        if(scale) res * N
        else res
      }
    )
  }

  if("norm" %in% distr) {
    norm.res = optim(par = c(hist.mean, sqrt(hist.var)),
                      fn = .hist.optim,
                      method = "L-BFGS-B",
                      .dist = "norm",
                      # fnscale = -1 does maximization and 1 for minimization
                      control = list(fnscale = 1),
                      lower = c(min(bin), 0.001),
                      upper = c(max(bin), (max(bin) - min(bin)) * 0.5))
    names(norm.res$par) = c("mean", "sd")
    norm.res$par <- as.list(norm.res$par)
    if(truncated) {
      norm.res$par$a =  min(bin) - 1e-10
      norm.res$par$b = max(bin) + 1e-10
    }
    norm.res$dist = "norm"
    norm.res$dens = function(x = NULL, scale = TRUE) {
      if(missing(x)) {
        x = bin
      }
      args = c(list(x = x), as.list(norm.res$par))
      res = do.call("dtnorm", args)
      if(scale) res * N
      else res
    }
    # Only
    if(norm.res$convergence == 0) rtn$norm = norm.res
    else warning("fit.distributions.optim: Normal distribution did not converge.")
  }


  if("gamma" %in% distr) {
    shape.init = hist.mean^2 / hist.var
    rate.init = hist.mean / hist.var
    gamma.res = optim(par = c(shape.init, rate.init),
                       fn = .hist.optim,
                       method = "L-BFGS-B",
                       .dist = "gamma",
                       control = list(fnscale = 1),
                       lower =c(0.001, 0.001))
    names(gamma.res$par) = c("shape", "rate")
    gamma.res$par <- as.list(gamma.res$par)
    if(truncated) {
      gamma.res$par$a =  min(bin) - 1e-10
      gamma.res$par$b = max(bin) + 1e-10
    }
    gamma.res$dist = "gamma"
    gamma.res$dens = function(x = NULL, scale = TRUE) {
      if(missing(x)) {
        x = bin
      }
      args = c(list(x = x), as.list(gamma.res$par))
      res = do.call("dtgamma", args)
      if(scale) res * N
      else res
    }
    if(gamma.res$convergence == 0) rtn$gamma = gamma.res
    else warning("fit.distributions.optim: Gamma distribution did not converge.")
  }

  rtn
}

# fit.histogram.metric.distributions = function(x, seg.start, seg.end, fit.mixtures = c("norm", "gamma"), optimal.uniform = TRUE)
