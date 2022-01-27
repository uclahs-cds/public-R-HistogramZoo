#' Convert a vector of points into a list of start/end points
#'
#' @param p a vector of indices
#'
#' @return A list with keys: start and end representing the indices
#' @export
#'
#' @examples
#' points.to.start.end(c(1,5,10))
points.to.start.end <- function(p) {
  if(length(p) <= 1) {
    stop("Need more than 1 point to compute start/end")
  }
  return.list = list(
    start = p[1:(length(p) - 1)]
  )
  return.list$end = c(p[2:(length(p) - 1)] - 1,  p[length(p)])
  return.list
}


#' Creates a new set of segments from a partition of points
#'
#' @param p a vector of points
#' @param max.gaps A list (or data frame) with start and end named columns
#' @param remove.short.segment
#'
#' @return
#' @export
#'
#' @examples
remove.max.gaps.agnostic = function(p, max.gaps, remove.short.segment = 1) {
  if(nrow(max.gaps) == 0) {
    # Return the original points as start/end data frame
    p.start.end = data.frame(
      start = p[1:(length(p) - 1)],
      end = p[2:length(p)])

    return(p.start.end)
  }
  p.seq = unlist(lapply(2:length(p), function(i) {
    seq(p[i - 1], p[i], by = 1)
    }))

  max.gaps.seq = unlist(lapply(1:nrow(max.gaps), function(i) {
    # TODO: Should this include the endpoints or no?
    seq(max.gaps[i, 1], max.gaps[i, 2], by = 1)
  }))

  p.no.maxgap = sort(unique(setdiff(p.seq, max.gaps.seq)))

  # https://stackoverflow.com/a/24837419
  # Take the segments and regroup into consecutive integers
  new.p.seq <- split(p.no.maxgap, cumsum(c(1, diff(p.no.maxgap) != 1)))
  new.p.seq <- new.p.seq[unlist(lapply(new.p.seq, length) > remove.short.segment)]

  # Only take the first and last consecutive numbers
  new.p <- lapply(new.p.seq, function(y) {
    list(start = y[1], end = y[length(y)])
  })

  do.call('rbind.data.frame', new.p)
}
