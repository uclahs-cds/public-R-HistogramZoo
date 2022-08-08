#' Creates a new set of segments from a partition of points
#'
#' @param start.end.points a data.frame with start and end named columns representing intervals
#' @param max.gaps a data.frame with start and end named columns representing gaps
#' @param remove.short.segment the minimum segment length to include in the final segments
#'
#' @return a data.frame of the intervals with max gaps removed
#' \describe{
#'     \item{start}{start indices of new intervals}
#'     \item{end}{end indices of new intervals}
#' }
#' @export
remove_max_gaps <- function(start.end.points, max.gaps, remove.short.segment = 0) {

  # Error checking
  if(!is.data.frame(start.end.points) | !all(c("start", "end") %in% colnames(start.end.points))){
    stop("start.end.points should be a data.frame with columns start and end")
  }
  if(!is_equal_integer(start.end.points[,"start"]) | !is_equal_integer(start.end.points[,"end"]) | !all(start.end.points[,"start"] <= start.end.points[,"end"])){
    stop("start.end.points start and end columns must represent functional intervals, i.e end >= start and end and start are numeric integers")
  }
  if(nrow(max.gaps) == 0) {
    # Return the original points as start/end data frame
    return(start.end.points)
  }
  if(!is.data.frame(max.gaps) | !all(c("start", "end") %in% colnames(max.gaps))){
    stop("max.gaps should be a data.frame with columns start and end")
  }
  if(!is_equal_integer(max.gaps[,"start"]) | !is_equal_integer(max.gaps[,"end"]) | !all(max.gaps[,"start"] <= max.gaps[,"end"])){
    stop("max.gaps start and end columns must represent functional intervals, i.e end >= start and end and start are numeric integers")
  }
  if(!is_equal_integer(remove.short.segment) | length(remove.short.segment) != 1){
    stop("remove.short.segment must be a numeric integer of length 1")
  }
  
  # For each segment, create a sequence of consecutive integers
  seg.seq.list <- mapply(seq.int, from = start.end.points$start, to = start.end.points$end, by = 1, SIMPLIFY = FALSE)

  max.gaps.seq <- unlist(lapply(1:nrow(max.gaps), function(i) {
    # TODO: Should this include the endpoints or no?
    seq(max.gaps[i, 1], max.gaps[i, 2], by = 1)
  }))

  seg.seq.no.maxgap <- lapply(seg.seq.list, function(p.seq) {
    sort(unique(setdiff(p.seq, max.gaps.seq)))
  })

  # https://stackoverflow.com/a/24837419
  # Take the segments and regroup into consecutive integers
  new.p.seq <- lapply(seg.seq.no.maxgap, function(p.no.maxgap) {
    consec.p <- split(p.no.maxgap, cumsum(c(1, diff(p.no.maxgap) != 1)))
    # Remove short segments
    # TODO: Should we do this later? Do we need to join any segments first?
    consec.p[unlist(lapply(consec.p, length) > remove.short.segment)]
  })

  # Remove one level of the lists for a flat structure
  new.p.seq.consec <- unname(unlist(new.p.seq, recursive = FALSE))

  # Extract the start and end points from the sequences
  new.p.start.end <- lapply(new.p.seq.consec, function(p.seq) {
    list('start' = p.seq[1], 'end' = p.seq[length(p.seq)])
  })

  res <- do.call('rbind.data.frame', new.p.start.end)
  rownames(res) <- NULL

  return(res)
}
