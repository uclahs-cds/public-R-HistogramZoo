
#' Fit the model parameters by optimizing a histogram metric
#'
#' @param freq TODO
#' @param metric TODO
#' @param truncated TODO
#' @param distr TODO
#' @export
#' @importFrom DEoptim DEoptim
fit_distributions = function(
  freq,
  metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
  truncated = FALSE,
  distr = c("norm", "gamma", "unif")) {

  # Matching arguments
  metric = match.arg(metric, several.ok = TRUE)
  distr = match.arg(distr, several.ok = TRUE)

  # Initializing Data
  L = length(freq)
  N = sum(freq)
  bin = 1:L
  # Used for initializing parameter estimations for norm and gamma
  # hist.mean = sum(freq * bin) / N
  # hist.var = sum(freq * (bin - hist.mean)^2) / N
  # shape.init = hist.mean^2 / hist.var
  # rate.init = hist.mean / hist.var

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
    if(is.na(res) || res == -Inf || res == Inf) {
      res = Inf
    }
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
        "par" = NULL,
        "dist" = "unif",
        "metric" = met,
        "value" = metric.func(freq, rep(unif.dens * N, L)),
        "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
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
      norm.optim = DEoptim::DEoptim(fn = .hist.optim,
                                    .dist = "norm",
                                    lower = c(min(bin), 0.001),
                                    upper = c(max(bin), (max(bin) - min(bin)) * 0.5),
                                    control = list(
                                    trace = FALSE, # Do not print results
                                    itermax = 500, # Iterations
                                    VTR = 10^-2 # At 1 %, stop optimizing
                                    ))
      names(norm.optim$optim$bestmem) = c("mean", "sd")
      norm.par = as.list(norm.optim$optim$bestmem)
      if(truncated) {
        norm.par$a =  min(bin) - 1e-10
        norm.par$b = max(bin) + 1e-10
      }
      rtn[[tag]] = list(
        "par" = norm.par,
        "dist" = "norm",
        "metric" = met,
        "value" = norm.optim$optim$bestval,
        "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
          if(missing(x)) {
            x = bin
          }
          args = c(list(x = x), as.list(mpar))
          res = do.call("dtnorm", args)
          if(scale) res * N
          else res
        }
      )
    }

    # Gamma Distribution
    if(distr == "gamma") {
      gamma.optim = DEoptim::DEoptim(fn = .hist.optim,
                                     .dist = "gamma",
                                     lower = c(0.001, 0.001),
                                     upper = c(L, L),
                                     control = list(
                                     trace = FALSE, # Do not print results
                                     itermax = 500, # Iterations
                                     VTR = 10^-2 # At 1 %, stop optimizing
                                     ))
      names(gamma.optim$optim$bestmem) = c("shape", "rate")
      gamma.par <- as.list(gamma.optim$optim$bestmem)
      if(truncated) {
        gamma.par$a =  min(bin) - 1e-10
        gamma.par$b = max(bin) + 1e-10
      }
      rtn[[tag]] = list(
        "par" = gamma.par,
        "dist" = "gamma",
        "metric" = met,
        "value" = gamma.optim$optim$bestval,
        "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
          if(missing(x)) {
            x = bin
          }
          args = c(list(x = x), as.list(mpar))
          res = do.call("dtgamma", args)
          if(scale) res * N
          else res
        }
      )
  }
  }

  rtn
}
