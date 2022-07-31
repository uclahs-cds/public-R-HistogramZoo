

#' create.residualplot creates a residual plot between fitted and observed data
#'
#' @param data 
#' @param distributions
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
create.residualplot = function(
  x,
  distributions,
  ...
){
  
  # Error checking
  
  # x
  if(!is.null(names(x))){
    labels.x = names(x)
    labels.x <- tryCatch(
      { labels.x = as.numeric(labels.x) },
      warning = function(cond) {
        message("Warning message:")
        message("Vector names are not coercible to numeric.")
        message("Using default indices.")
        return(1:length(x))
      }
    )
  } else {
    labels.x = 1:length(x)
  }
  plotting.data = data.frame("density" = x, "labels.x" = labels.x)
  
  # distribution
  # Extract fitted distributions of majority vote model
  if(!is.null(distributions)){
    models = distributions[['models']]
    majority_vote = lapply(models, `[[`, "majority.vote")
    distribution_plotting_data = lapply(majority_vote, function(m) {
      x = seq(m$seg.start, m$seg.end, by = 1)
      dens = m$dens(x = seq_along(x), mpar = m$par)
      data.frame("fitted" = dens, "labels.x" = labels.x[x])
    })
    distribution_plotting_data = do.call('rbind.data.frame', distribution_plotting_data)
    
  }
  
  # Calculating residuals
  plotting.data = merge(plotting.data, distribution_plotting_data, by = "labels.x", all = T)
  plotting.data$Residuals = plotting.data$density - plotting.data$fitted
  
  # Adding lines to help with visualization (Can potentially remove)
  plotting.chgpts = which(abs(diff(sign(plotting.data$Residuals))) == 2)
  
  # Plotting
  plt =  BoutrosLab.plotting.general::create.scatterplot(
    Residuals ~ labels.x,
    plotting.data,
    # Colour
    col = "black",
    # Lines & PCH
    type = c('p'),
    # Adding lines at changepoints
    abline.v = plotting.chgpts,
    abline.h = 0,
    abline.lty = "dotted",
    abline.col = "lightgrey",
    abline.lwd = 0.01
  )
  return(plt)
  
  # Return plt
  return(plt)
}