

#' create_trackplot creates a horizontal plot with rows representing separate tracks
#'
#' @param track_data A dataframe detailing segments
#' \describe{
#'     \item{start}{the start position of the segment, based on histogram vector index}
#'     \item{end}{the end position of the segment, based on histogram vector index}
#'     \item{`metric_id`}{values; factor or continuous data for plotting}
#'     \item{`row_id`}{character or factor; identifier for the row}
#'     \item{colour}{colour for the row}
#' }
#'
#' TODO: Add attributes of x for labels and stuff?
#'
#' @return
#' @export
#'
#' @examples
create_trackplot = function(
  track_data
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
