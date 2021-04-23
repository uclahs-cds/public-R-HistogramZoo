
extract.distribution.parameters = function(
  mod,
  x,
  scalefactor
){
  
  extract.distribution.parameters.helper = function(m){
    
    if(class(m) == "mixEM"){
      mixfit.params <- length(m$mu)*3 # for lambda, mu, and sigma params
      
      # For Residuals
      bin.data <- table(x)
      dens <- dnorm_mixture(as.numeric(names(bin.data)), m)
      dens.scale <- dens * scalefactor
      fit.residuals <- (dens.scale - as.integer(bin.data))^2
      
      results <- data.frame(
        "dist" = "norm_mixture",
        "loglikelihood" = m$loglik,
        "aic" = -2 * m$loglik + 2 * mixfit.params,
        "bic" = -2 * m$loglik + mixfit.params * log(length(x)),
        "mse" = mean(fit.residuals),
        "params" = dput.str(m[c("mu", "sigma", "lambda")]),
        stringsAsFactors = F
      )
    } else {
      
      # Residuals
      bin.data <- table(x)
      params <- c(as.list(m$estimate), as.list(m$fix.arg))
      distname <- m$distname
      ddistname <- paste0("d", distname)
      call.params <- c(list(x = as.numeric(names(bin.data))), as.list(params))
      dens <- do.call(ddistname, call.params)
      dens.scale <- dens * scalefactor
      
      fit.residuals <- (dens.scale - as.integer(bin.data))^2
      
      data.frame(
        "dist" = summary(m)$distname,
        "loglikelihood" = summary(m)$loglik,
        "aic" = summary(m)$aic,
        "bic" = summary(m)$bic,
        "mse" = mean(fit.residuals),
        # Text representation of the parameters
        params = dput.str(params),
        stringsAsFactors = F)
    }
  }
  
  do.call(rbind, lapply(mod, extract.distribution.parameters.helper))
}