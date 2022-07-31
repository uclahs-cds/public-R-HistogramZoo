

#' create_residualplot creates a residual plot between fitted and observed data
#'
#' @param histogram_obj a `Histogram` or `HistogramFit` object
#' @param model_name One of the metrics used to fit models (e.g. Jaccard) and "consensus", default consensus
#' @param ... Additional parameters for the BoutrosLab.plotting.general function create.scatterplot
#'
#' @return Residual scatterplot, a Trellis object. For further details, see the 'Lattice' R package.
#' @export
#'
#' @examples
create_residualplot <- function(histogram_obj, model_name){
  UseMethod('create_residualplot')
}

create_residualplot.HistogramFit = function(
  histogram_obj,
  model_name,
  ...
){

  # Error checking
  stopifnot(inherits(histogram_obj, "HistogramFit"))
  model_name <- match.arg(model_name, c("consensus", histogram_obj$histogram_metric))

  # Extracting histogram_data
  x <- histogram_obj$histogram_data
  xaxis.labels <- generate_interval_labels(x$interval_start, x$interval_end)
  plotting.data <- data.frame(
    "density" = x,
    "labels.x" = 1:length(x)
  )

  # distribution
  models <- distributions[['models']]
  mods <- lapply(models, `[[`, model_name)
  distribution_plotting_data <- lapply(mods, function(m) {
    x <- seq(m$seg.start, m$seg.end, by = 1)
    dens <- m$dens(x = seq_along(x), mpar = m$par)
    return(
      data.frame("fitted" = dens, "labels.x" = x)
    )
  })
  distribution_plotting_data <- do.call('rbind.data.frame', distribution_plotting_data)

  # Calculating residuals
  plotting.data <- merge(plotting.data, distribution_plotting_data, by = "labels.x", all = T)
  plotting.data$Residuals <- plotting.data$density - plotting.data$fitted

  # Adding lines to help with visualization (Can potentially remove)
  plotting.chgpts <- which(abs(diff(sign(plotting.data$Residuals))) == 2)

  # Plotting
  plt <-  BoutrosLab.plotting.general::create.scatterplot(
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
    abline.lwd = 0.01,
    ...
  )

  # Return plt
  return(plt)
}
