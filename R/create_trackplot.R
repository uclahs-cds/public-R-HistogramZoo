

#' create.trackplot creates a horizontal plot with rows representing separate tracks
#'
#' @param x A vector of histogram coverage
#' @param track_data A dataframe detailing segments
#' \describe{
#'     \item{start}{the start position of the segment, based on histogram vector index}
#'     \item{end}{the end position of the segment, based on histogram vector index}
#'     \item{`metric_id`}{values; factor or continuous data for plotting}
#'     \item{`row_id`}{character or factor; identifier for the row}
#' }
#' @param metric_id The name of the column containing the data used for plotting
#' @param row_id The name of the column containing the data used to indicate the row of the track
#' @param colour.scheme A colour scheme for the data
#' @param ...
#'
#' TODO: Add attributes of x for labels and stuff?
#'
#' @return
#' @export
#'
#' @examples
create.trackplot = function(
  x = NULL,
  track_data,
  metric_id,
  row_id,
  colour.scheme = distribution_colours,
  scale.colours = c("continuous", "discrete"),
  ...
){

  # Error checking
  scale.colours = match.arg(scale.colours, several.ok = FALSE)
  if(scale.colours == "continuous"){
    # Check that this is a colour vector of length 2
    # Otherwise use the first two colours
  } else if (scale.colours == "discrete"){
    # Check that there are enough colours for each of the categories
  }

  # Rows
  rows = unique(track_data[,row_id])
  nrows = length(rows)

  # Columns
  if(!is.null(x)){ # Add xlimits, figure out how to do this better
    x.limits = c(1, length(x))
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
  }

  # Initializing matrix
  initial.matrix = matrix(
    NA,
    ncol = nrows,
    nrow = x.limits[2],
    dimnames = list(labels.x, rows)
  )

  # Encoding for discrete colours
  if(scale == "discrete"){
    track_data[,metric_id] = factor(track_data[,metric_id], levels = names(colour.scheme))
  }

  # Filling in the matrix
  for(i in 1:nrow(track_data) ){
    initial.matrix[track_data[i, "start"]:track_data[i, "end"], track_data[i, row_id]] <- track_data[i, metric_id]
  }

  # A hack for when there's only one row
  if(nrows == 1){
    initial.matrix = cbind(initial.matrix, initial.matrix)
  }

  # Colour scale # Add this as a default to a user-set parameter
  colour.at = if(scale == "discrete") seq(0, length(colour.scheme), 1) + 0.5 else seq(min(track_data[, metric_id]), max(track_data[, metric_id]), length.out = 10)

  # Making a plot
  plt = BoutrosLab.plotting.general::create.heatmap(
    initial.matrix,
    clustering.method = 'none',
    # Colours
    colour.scheme = colour.scheme,
    at = colour.at,
    fill.colour = "white",
    # XX
    print.colour.key = F
  )

  # Return plt
  return(plt)
}
