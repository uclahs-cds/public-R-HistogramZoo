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
