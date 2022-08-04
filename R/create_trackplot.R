

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
#' @param row_id character, column name used to indicate the track row
#' @param metric_id character, column name used to scale colour
#' @param xlimits numeric vector of length 2 indicating the start and end indices of the track respectively
#' @inheritParams BoutrosLab.plotting.general::create.heatmap
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
create_trackplot = function(
  track_data,
  row_id,
  metric_id,
  xlimits = c(min(track_data$start), max(track_data$end)),
  # BPG params for colours
  colour.scheme = c("slateblue4", "violetred3"),
  total.colours = 99,
  colour.centering.value = 0,
  colour.alpha = 1,
  fill.colour = 'white',
  at = NULL,
  print.colour.key = TRUE,
  colourkey.labels.at = NULL,
  colourkey.labels = NULL,
  # BPG defaults
  filename = NULL,
  main = list(label = ''),
  main.just = "center",
  main.x = 0.5,
  main.y = 0.5,
  main.cex = 3,
  yaxis.lab = NULL,
  xaxis.lab = NULL,
  xaxis.lab.top = NULL,
  xaxis.cex = 1.5,
  xaxis.top.cex = NULL,
  yaxis.cex = 1.5,
  xlab.cex = 2,
  ylab.cex = 2,
  xlab.top.label = NULL,
  xlab.top.cex = 2,
  xlab.top.col = 'black',
  xlab.top.just = "center",
  xlab.top.x = 0.5,
  xlab.top.y = 0,
  xat = TRUE,
  xat.top = NULL,
  yat = TRUE,
  xaxis.tck = NULL,
  xaxis.top.tck = NULL,
  yaxis.tck = NULL,
  xaxis.col = 'black',
  yaxis.col = 'black',
  col.pos = NULL,
  row.pos = NULL,
  cell.text = '',
  text.fontface = 1,
  text.cex = 1,
  text.col = 'black',
  text.position = NULL,
  text.offset = 0,
  text.use.grid.coordinates = TRUE,
  colourkey.cex = 3.6,
  xaxis.rot = 90,
  xaxis.rot.top = 90,
  yaxis.rot = 0,
  xlab.label = '' ,
  ylab.label = '',
  xlab.col = 'black',
  ylab.col = 'black',
  axes.lwd = 2,
  gridline.order = 'h',
  grid.row = FALSE,
  grid.col = FALSE,
  force.grid.row = FALSE,
  force.grid.col = FALSE,
  grid.limit = 50,
  row.lines = seq(0, ncol(x), 1) + 0.5,
  col.lines = seq(0, nrow(x), 1) + 0.5,
  top.padding = 0.1,
  bottom.padding = 0.5,
  right.padding = 0.5,
  left.padding = 0.5,
  x.alternating = 1,
  shrink = 1,
  row.colour = 'black',
  col.colour = 'black',
  row.lwd = 1,
  col.lwd = 1,
  grid.colour = NULL,
  grid.lwd = NULL,
  width = 6,
  height = 6,
  size.units = 'in',
  resolution = 1600,
  enable.warnings = FALSE,
  xaxis.fontface = 'bold',
  yaxis.fontface = 'bold',
  input.colours = FALSE,
  axis.xlab.padding = 0.1,
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
    stop("row_id needs to be factor, character or numeric")
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
  
  # Rows
  rows = sort(unique(track_data[,row_id]))
  nrows = length(rows)
  
  # Columns
  cols = seq(xlimits[1], xlimits[2])
  ncols = length(cols)

  # Initializing matrix
  initial.matrix = matrix(
    NA,
    nrow = nrows,
    ncol = ncols,
    dimnames = list(rows, cols)
  )

  # Filling in the matrix
  for(i in 1:nrow(track_data)){
    itv_start = max(c(1, track_data[i, "start"] - xlimits[1] + 1))
    itv_end = min(c(xlimits[2] - xlimits[1] + 1, track_data[i, "end"] - xlimits[1] + 1))
    initial.matrix[track_data[i, row_id], itv_start:itv_end] <- track_data[i, metric_id]
  }

  # Making a plot
  plt = BoutrosLab.plotting.general::create.heatmap(
    initial.matrix,
    clustering.method = 'none',
    same.as.matrix = T,
    # Colours
    colour.scheme = colour.scheme,
    total.colours = total.colours,
    colour.centering.value = colour.centering.value,
    colour.alpha = colour.alpha,
    fill.colour = fill.colour,
    at = at,
    print.colour.key = print.colour.key,
    colourkey.labels.at = colourkey.labels.at,
    colourkey.labels = colourkey.labels,
    # BPG defaults
    filename = filename,
    main = main,
    main.just = main.just,
    main.x = main.x,
    main.y = main.y,
    main.cex = main.cex,
    yaxis.lab = yaxis.lab,
    xaxis.lab = xaxis.lab,
    xaxis.lab.top = xaxis.lab.top,
    xaxis.cex = xaxis.cex,
    xaxis.top.cex = xaxis.top.cex,
    yaxis.cex = yaxis.cex,
    xlab.cex = xlab.cex,
    ylab.cex = ylab.cex,
    xlab.top.label = xlab.top.label,
    xlab.top.cex = xlab.top.cex,
    xlab.top.col = xlab.top.col,
    xlab.top.just = xlab.top.just,
    xlab.top.x = xlab.top.x,
    xlab.top.y = xlab.top.y,
    xat = xat,
    xat.top = xat.top,
    yat = yat,
    xaxis.tck = xaxis.tck,
    xaxis.top.tck = xaxis.top.tck,
    yaxis.tck = yaxis.tck,
    xaxis.col = xaxis.col,
    yaxis.col = yaxis.col,
    col.pos = col.pos,
    row.pos = row.pos,
    cell.text = cell.text,
    text.fontface = text.fontface,
    text.cex = text.cex,
    text.col = text.col,
    text.position = text.position,
    text.offset = text.offset,
    text.use.grid.coordinates = text.use.grid.coordinates,
    colourkey.cex = colourkey.cex,
    xaxis.rot = xaxis.rot,
    xaxis.rot.top = xaxis.rot.top,
    yaxis.rot = yaxis.rot,
    xlab.label = xlab.label,
    ylab.label = ylab.label,
    xlab.col = xlab.col,
    ylab.col = ylab.col,
    axes.lwd = axes.lwd,
    gridline.order = gridline.order,
    grid.row = grid.row,
    grid.col = grid.col,
    force.grid.row = force.grid.row,
    force.grid.col = force.grid.col,
    grid.limit = grid.limit,
    row.lines = row.lines,
    col.lines = col.lines,
    top.padding = top.padding,
    bottom.padding = bottom.padding,
    right.padding = right.padding,
    left.padding = left.padding,
    x.alternating = x.alternating,
    shrink = shrink,
    row.colour = row.colour,
    col.colour = col.colour,
    row.lwd = row.lwd,
    col.lwd = col.lwd,
    grid.colour = grid.colour,
    grid.lwd = grid.lwd,
    width = width,
    height = height,
    size.units = size.units,
    resolution = resolution,
    enable.warnings = enable.warnings,
    xaxis.fontface = xaxis.fontface,
    yaxis.fontface = yaxis.fontface,
    input.colours = input.colours,
    axis.xlab.padding = axis.xlab.padding,
    ...
  )

  # Return plt
  return(plt)
}
