extract.distribution.parameters = function(
  mod,
  x
){

  extract.distribution.parameters.helper = function(m){

    if(class(m) == "mixEM"){

      # Extracting Model Parameters
      mixfit.params <- length(m$mu)*3 # for lambda, mu, and sigma params

      # Density of Observed Data
      bin.data <- table(x)
      scaled.observations = as.integer(bin.data)/sum(as.integer(bin.data))

      # Calculating Density
      dens <- dnorm_mixture(as.numeric(names(bin.data)), m)

      # Calculating Residuals
      fit.residuals = (dens - scaled.observations)^2

      # Jaccard
      overlap = pmin(a = scaled.observations, b = dens, na.rm = T)
      union = pmax(a = scaled.observations, b = dens, na.rm = T)
      jc = sum(overlap)/sum(union)

      results <- data.frame(
        "dist" = "norm_mixture",
        "loglikelihood" = m$loglik,
        "aic" = -2 * m$loglik + 2 * mixfit.params,
        "bic" = -2 * m$loglik + mixfit.params * log(length(x)),
        "mse" = mean(fit.residuals),
        "jc" = jc,
        "params" = dput.str(m[c("mu", "sigma", "lambda")]),
        stringsAsFactors = F
      )
    } else {

      # Density of Observed Data
      bin.data <- table(x)
      scaled.observations = as.integer(bin.data)/sum(as.integer(bin.data))

      # Extracting Model Parameters & Calculating Density
      params <- c(as.list(m$estimate), as.list(m$fix.arg))
      distname <- m$distname
      ddistname <- paste0("d", distname)
      call.params <- c(list(x = as.numeric(names(bin.data))), as.list(params))
      dens <- do.call(ddistname, call.params)

      # Calculating Residuals
      fit.residuals = (dens - scaled.observations)^2

      # Jaccard
      overlap = pmin(a = scaled.observations, b = dens, na.rm = T)
      union = pmax(a = scaled.observations, b = dens, na.rm = T)
      jc = sum(overlap)/sum(union)

      data.frame(
        "dist" = summary(m)$distname,
        "loglikelihood" = summary(m)$loglik,
        "aic" = summary(m)$aic,
        "bic" = summary(m)$bic,
        "mse" = mean(fit.residuals),
        "jc" = jc,
        params = dput.str(params),
        stringsAsFactors = F)
    }
  }

  do.call(rbind, lapply(mod, extract.distribution.parameters.helper))
}
