
# Defining plotting parameters
distribution_colours = c(
  "coverage" = "black",
  "norm" = "darkorange",
  "gamma" = "chartreuse4",
  "unif" = "darkorchid4"
)

distribution_lwd = c(
  "coverage" = 1,
  "norm" = 2.5,
  "gamma" = 2.5,
  "unif" = 2.5
)

#' create_coverage_plot creates a histogram/coverage plot with the potential to add annotations
#'
#' @param x
#' @param points
#' @param distributions
#' @param fitted_data
#'
#' @return
#' @export
#'
#' @examples
create.coverageplot = function(
  x,
  points = NULL,
  distributions = NULL,
  fitted_data = NULL,
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
  plotting.data = data.frame(x, labels.x)

  # points
  if(!is.null(points)){
    points.x = labels.x[points]
    points.y = x[points]
  } else {
    points.x = NULL
    points.y = NULL
  }

  # distribution
  if(!is.null(distributions)){
    # TODO
    # plotting.data
    # distribution_groups
    # colour_vector
    # lwd_vector
  } else {
    distribution_groups = rep(1, length(x))
    color_vector = "black"
    lwd_vector = 1
  }
  
  # Plotting
  plt = BoutrosLab.plotting.general::create.scatterplot(
    x ~ labels.x,
    data = plotting.data,
    # Groups
    groups = distribution_groups,
    col = color_vector,
    lwd = lwd_vector,
    # Lines & PCH
    type = c('a'),
    # Adding extra points
    add.points = T,
    points.x = points.x,
    points.y = points.y,
    points.pch = 19,
    points.col = 'red',
    ...
  )
  
  # Returning plotting object
  return(plt)
}