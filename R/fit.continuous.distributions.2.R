
#' Fit the model parameters by optimizing a histogram metric
fit.distributions.optim = function(freq, metric = c("jaccard", "intersection", "ks", "mse", "chisq"), truncated = FALSE, distr = c("norm", "gamma", "unif")) {

  # Matching arguments
  metric = match.arg(metric, several.ok = TRUE)
  distr = match.arg(distr, several.ok = TRUE)

  # Initializing Data
  L = length(freq)
  N = sum(freq)
  bin = 1:length(freq)
  hist.mean = sum(freq * bin) / N
  hist.var = sum(freq * (bin - hist.mean)^2) / N
  shape.init = hist.mean^2 / hist.var
  rate.init = hist.mean / hist.var

  # Optimization Function
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

  # Initializing Todo & Results
  todo = expand.grid(metric, distr, stringsAsFactors = F)
  rtn = list()
  for(i in 1:nrow(todo)){

    # Extracting Metric & Distribution
    met = todo[i,1]
    distr = todo[i,2]
    tag = paste0(met, ".", distr)

    # Get one of the metrics from histogram.distances
    metric.func = get(paste('histogram', met, sep = "."))

    # Uniform Distribution
    if(distr == "unif") {
      # Add uniform distribution
      unif.dens = 1 / (max(bin) - min(bin))
      rtn[[tag]] = list(
        "dist" = "unif",
        "metric" = met,
        "value" = metric.func(freq, rep(unif.dens * N, L)),
        "dens" = function(x = NULL, scale = TRUE) {
          if(missing(x)) {
            x = bin
          }
          res = ifelse(x >= min(bin) & x <= max(bin), unif.dens, 0)
          if(scale) res * N
          else res
        }
      )
    }
    # Normal Distribution
    if(distr == "norm") {
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
      norm.res$metric = met
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
      if(norm.res$convergence == 0) rtn[[tag]] = norm.res
      else warning("fit.distributions.optim: Normal distribution did not converge.")
    }
    # Gamma Distribution
    if(distr == "gamma") {
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
      gamma.res$metric = met
      gamma.res$dens = function(x = NULL, scale = TRUE) {
        if(missing(x)) {
          x = bin
        }
        args = c(list(x = x), as.list(gamma.res$par))
        res = do.call("dtgamma", args)
        if(scale) res * N
        else res
      }
      if(gamma.res$convergence == 0) rtn[[tag]] = gamma.res
      else warning("fit.distributions.optim: Gamma distribution did not converge.")
    }
  }

  rtn
}
