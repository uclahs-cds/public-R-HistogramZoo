
#' Stacks together multiple plots
#'
#' @inheritParams BoutrosLab.plotting.general::create.multipanelplot 
#' 
#' @return a Trellis object. For further details, see the 'Lattice' R package.
#' @export 
#'
#' @examples \dontrun{
#' x = rnorm(10000, mean = 100, sd = 50)
#' x = observations_to_histogram(round(x), histogram_bin_width = 5)
#' results = segment_and_fit(x, eps = 1)
#' cvg_plt = create_coverageplot(
#'   results
#' )
#' resid_plt = create_residualplot(
#' results
#' )
#' 
#' create_layerplot(
#' plot.objects = list(cvg_plt, resid_plt)
#' )
#' }
create_layerplot <- function(
    plot.objects = NULL,
    filename = NULL,
    height = 10,
    width = 10,
    resolution = 1000,
    plot.objects.heights = rep(1, length(plot.objects)),
    main = '',
    main.x = 0.5,
    main.y = 0.5,
    x.spacing = 0, 
    y.spacing = 0,
    xlab.label = '',
    xlab.cex = 2,
    ylab.label = '',
    ylab.label.right = '',
    ylab.cex = 2,
    main.cex = 3,
    legend = NULL,
    left.padding = 0,
    ylab.axis.padding = 0,
    xlab.axis.padding = 0,
    bottom.padding = 0,
    top.padding = 0,
    right.padding = 0,
    left.legend.padding = 2,
    right.legend.padding = 2, 
    bottom.legend.padding = 2, 
    top.legend.padding = 2,
    size.units = 'in',
    enable.warnings = FALSE
){
  
  plt <- BoutrosLab.plotting.general::create.multipanelplot(
    plot.objects = plot.objects,
    filename = filename,
    height = height,
    width = width,
    resolution = resolution,
    plot.objects.heights = plot.objects.heights,
    plot.objects.widths = 1,
    layout.width = 1,
    layout.height = length(plot.objects),
    main = main,
    main.x = main.x,
    main.y = main.y,
    x.spacing = x.spacing, 
    y.spacing = y.spacing,
    xlab.label = xlab.label,
    xlab.cex = xlab.cex,
    ylab.label = ylab.label,
    ylab.label.right = ylab.label.right,
    ylab.cex = ylab.cex,
    main.cex = main.cex,
    legend = legend,
    left.padding = left.padding,
    ylab.axis.padding = ylab.axis.padding,
    xlab.axis.padding = xlab.axis.padding,
    bottom.padding = bottom.padding,
    top.padding = top.padding,
    right.padding = right.padding,
    layout.skip = rep(FALSE, length(plot.objects)),
    left.legend.padding = left.legend.padding,
    right.legend.padding = right.legend.padding, 
    bottom.legend.padding = bottom.legend.padding, 
    top.legend.padding = top.legend.padding,
    size.units = size.units,
    enable.warnings = enable.warnings
  )
  
  return(plt)
}