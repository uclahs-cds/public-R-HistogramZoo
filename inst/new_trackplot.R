# https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
#' Verify if a vector contains R colours
#' @param x A vector of potential colurs
#' @return logical vector indicating which elements are colours
are_colours <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(grDevices::col2rgb(X)), 
             error = function(e) FALSE)
  })
}

#' Generates row names from row data
#' @param row_data a vector of row data
#' @return a vector of row names
generate_row_ids <- function(row_data){
  if(is.factor(row_data)){
    row_names <- levels(row_data)
  } else if (is.character(row_data)){
    row_names <- unique(row_data)
  } else {
    row_names <- sort(unique(row_data))
  }
  return(row_names)
}

#' create_trackplot creates a horizontal plot with rows representing separate tracks
#'
#' @param track_data a dataframe detailing segments
#' \describe{
#'     \item{start}{the start position of the segment, based on histogram vector index}
#'     \item{end}{the end position of the segment, based on histogram vector index}
#'     \item{`metric_id`}{values; numeric data for plotting}
#'     \item{`row_id`}{character or factor; identifier for the row}
#'     \item{`colour_id`}{colour for the row}
#' }
#' @param row_id character, column name used to indicate the track row
#' @param metric_id character, column name of metric column
#' @param colour_id character, optional, column name for the colour of the track
#' @param colour_scheme a two colour scheme to scale numeric data in the `metric_id` column, default white, black
#' @param alpha alpha (transparency) of tracks, default 1
#' @param xlimits numeric vector of length 2 indicating the start and end indices of the track respectively
#' @inheritParams BoutrosLab.plotting.general::create.scatterplot
#' @param ... Additional arguments to be passed to BoutrosLab.plotting.general::create.scatterplot
#'
#' @return Track plot, a Trellis object. For further details, see the 'Lattice' R package.
#'
#' @export
#'
#' @examples \dontrun{
#' x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
#' results = segment_and_fit(x, eps = 0.005)
#' results_table = summarize_results(results)
#' create_trackplot(
#' track_data = results_table,
#' row_id = "peak_id",
#' metric_id = "value"
#' )
#' }
create_trackplot <- function(
    track_data,
    row_id,
    metric_id,
    colour_id = NULL,
    colour_scheme = c("red", "blue"),
    alpha = 1,
    xlimits = c(min(track_data$start), max(track_data$end)),
    # BPG defaults
    filename = NULL,
    main = NULL,
    main.just = 'center',
    main.x = 0.5,
    main.y = 0.5,
    main.cex = 3,
    xlab.label = "Interval Coordinates",
    ylab.label = "Row",
    xlab.cex = 2,
    ylab.cex = 2,
    xlab.col = 'black',
    ylab.col = 'black',
    xlab.top.label = NULL,
    xlab.top.cex = 2,
    xlab.top.col = 'black',
    xlab.top.just = 'center',
    xlab.top.x = 0.5,
    xlab.top.y = 0,
    xat = TRUE,
    xaxis.lab = NA,
    yaxis.lab = generate_row_ids(track_data[,row_id]),
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xaxis.rot = 0,
    yaxis.rot = 0,
    xaxis.fontface = 'bold',
    yaxis.fontface = 'bold',
    xaxis.col = 'black',
    yaxis.col = 'black',
    xaxis.tck = c(1,1),
    yaxis.tck = c(1,1),
    add.grid = FALSE,
    xgrid.at = xat,
    ygrid.at = yat,
    grid.colour = NULL,
    axes.lwd = 1,
    key = list(text = list(lab = c(''))),
    legend = NULL,
    top.padding = 0.1,
    bottom.padding = 0.7,
    right.padding = 0.1,
    left.padding = 0.5,
    key.top = 0.1,
    key.left.padding = 0,
    ylab.axis.padding = 1,
    axis.key.padding = 1,
    add.axes = FALSE,
    axes.lty = 'dashed',
    abline.h = NULL,
    abline.v = NULL,
    abline.col = 'black',
    abline.lwd = 1,
    abline.lty = 1,
    add.points = FALSE,
    points.x = NULL,
    points.y = NULL,
    points.pch = 19,
    points.col = 'black',
    points.col.border = 'black',
    points.cex = 1,
    add.line.segments = FALSE,
    line.start = NULL,
    line.end = NULL,
    line.col = 'black',
    line.lwd = 1,
    height = 6,
    width = 6,
    size.units = 'in',
    resolution = 1600,
    enable.warnings = FALSE,
    ...
){
  
  # Error checking
  stopifnot(is.data.frame(track_data))
  if(!row_id %in% colnames(track_data)){
    stop("row_id not found in track_data columns")
  }
  if(!metric_id %in% colnames(track_data)){
    stop("metric_id not found in track_data columns")
  }
  if(!is.numeric(track_data[,metric_id])){
    stop("metric_id column needs to be numeric")
  }
  if( !is.factor(track_data[,row_id]) & !is.character(track_data[,row_id]) & !is.numeric(track_data[,row_id]) ){
    stop("row_id needs to be factor, character or integer")
  }
  if(!is_equal_integer(xlimits)){
    stop("xlimits must be integer")
  }
  if(!is_equal_integer(track_data[,"start"]) | ! is_equal_integer(track_data[,"end"]) ){
    stop("start and end must be integer")
  }
  if(!all(track_data[,"start"] <= track_data[,"end"])){
    stop("start and end must represent valid intervals")
  }
  if(!all(track_data[,"start"] >= xlimits[1]) | !all(track_data[,"end"] <= xlimits[2])){
    warning("some intervals are truncated by xlimits")
  }
  # Colour scheme
  if(!is.null(colour_id) & !all(are_colours(track_data[,colour_id]))){
    stop("all values in `colour_id` must represent real colours")
  }
  
  # Row
  row_names <- generate_row_ids(track_data[,row_id])
  if(is.factor(track_data[,row_id])){
    rows <- as.numeric(track_data[,row_id])
  } else if (is.character(track_data[,row_id])){
    rows <- factor(track_data[,row_id], levels = row_names)
    rows <- as.numeric(rows)
  } else {
    rows <- rank(track_data[,row_id], ties.method = "min")
  }
  nrows <- length(unique(rows))
  
  # Colour scheme
  if(!is.null(colour_id)){
    colour_vec <- track_data[,colour_id]
  } else {
    # TODO: potentially make this more complicated to involve better scaling, etc. 
    scaled_metric <- (track_data[,metric_id] - min(track_data[,metric_id]))/(max(track_data[,metric_id]) - min(track_data[,metric_id]))
    colour_ramp <- grDevices::colorRamp(colour_scheme)
    colour_vec <- grDevices::rgb(colour_ramp(scaled_metric), maxColorValue = 255)
  }
  
  # Generating plot
  plt <- BoutrosLab.plotting.general::create.scatterplot(
    # Dummy data
    formula = start ~ end,
    data = track_data,
    col = "transparent",
    # Limits
    xlimits = xlimits,
    ylimits = c(0, nrows) + 0.5,
    yat = seq(1, nrows, 1),
    # Bread and butter
    add.rectangle = T,
    xleft.rectangle = track_data$start,
    ybottom.rectangle = rows - 0.5,
    xright.rectangle = track_data$end,
    ytop.rectangle = rows + 0.5,
    # Colour scheme
    col.rectangle = colour_vec,
    alpha.rectangle = alpha,
    # BPG Defaults
    filename = filename,
    main = main,
    main.just = main.just,
    main.x = main.x,
    main.y = main.y,
    main.cex = main.cex,
    xlab.label = xlab.label,
    ylab.label = ylab.label,
    xlab.cex = xlab.cex,
    ylab.cex = ylab.cex,
    xlab.col = xlab.col,
    ylab.col = ylab.col,
    xlab.top.label = xlab.top.label,
    xlab.top.cex = xlab.top.cex,
    xlab.top.col = xlab.top.col,
    xlab.top.just = xlab.top.just,
    xlab.top.x = xlab.top.x,
    xlab.top.y = xlab.top.y,
    xat = xat,
    xaxis.lab = xaxis.lab,
    yaxis.lab = yaxis.lab,
    xaxis.cex = xaxis.cex,
    yaxis.cex = yaxis.cex,
    xaxis.rot = xaxis.rot,
    yaxis.rot = yaxis.rot,
    xaxis.fontface = xaxis.fontface,
    yaxis.fontface = xaxis.fontface,
    xaxis.col = xaxis.col,
    yaxis.col = yaxis.col,
    xaxis.tck = xaxis.tck,
    yaxis.tck = yaxis.tck,
    add.grid = add.grid,
    xgrid.at = xgrid.at,
    ygrid.at = ygrid.at,
    grid.colour = grid.colour,
    axes.lwd = axes.lwd,
    key = key,
    legend = legend,
    top.padding = top.padding,
    bottom.padding = bottom.padding,
    right.padding = right.padding,
    left.padding = left.padding,
    key.top = key.top,
    key.left.padding = key.left.padding,
    ylab.axis.padding = ylab.axis.padding,
    axis.key.padding = axis.key.padding,
    add.axes = add.axes,
    axes.lty = axes.lty,
    abline.h = abline.h,
    abline.v = abline.v,
    abline.col = abline.col,
    abline.lwd = abline.lwd,
    abline.lty = abline.lty,
    add.points = add.points,
    points.x = points.x,
    points.y = points.y,
    points.pch = points.pch,
    points.col = points.col,
    points.col.border = points.col.border,
    points.cex = points.cex,
    add.line.segments = add.line.segments,
    line.start = line.start,
    line.end = line.end,
    line.col = line.col,
    line.lwd = line.lwd,
    height = height,
    width = width,
    size.units = size.units,
    resolution = resolution,
    enable.warnings = enable.warnings,
    ...
  )
  
  # Return plt
  return(plt)
}
