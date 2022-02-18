#' Creates a new set of segments from a partition of points
#'
#' @param p a vector of points
#' @param max.gaps A list (or data frame) with start and end named columns
#' @param remove.short.segment The minimum segment length to include in the final segments
#'
#' @return TODO
#' @export
remove.max.gaps.agnostic = function(p, max.gaps, remove.short.segment = 0) {
  # Maybe move this out of the function?
  start.end.points = index.to.start.end(p)
  if(nrow(max.gaps) == 0) {
    # Return the original points as start/end data frame
    return(start.end.points)
  }

  # For each segment, create a sequence of consecutive integers
  seg.seq.list = mapply(seq.int, from = start.end.points$start, to = start.end.points$end, by = 1, SIMPLIFY = FALSE)

  max.gaps.seq = unlist(lapply(1:nrow(max.gaps), function(i) {
    # TODO: Should this include the endpoints or no?
    seq(max.gaps[i, 1], max.gaps[i, 2], by = 1)
  }))

  seg.seq.no.maxgap = lapply(seg.seq.list, function(p.seq) {
    sort(unique(setdiff(p.seq, max.gaps.seq)))
  })

  # https://stackoverflow.com/a/24837419
  # Take the segments and regroup into consecutive integers
  new.p.seq = lapply(seg.seq.no.maxgap, function(p.no.maxgap) {
    consec.p = split(p.no.maxgap, cumsum(c(1, diff(p.no.maxgap) != 1)))
    # Remove short segments
    # TODO: Should we do this later? Do we need to join any segments first?
    consec.p[unlist(lapply(consec.p, length) > remove.short.segment)]
  })

  # Remove one level of the lists for a flat structure
  new.p.seq.consec = unname(unlist(new.p.seq, recursive = FALSE))

  # Extract the start and end points from the sequences
  new.p.start.end = lapply(new.p.seq.consec, function(p.seq) {
    list(start = p.seq[1], end = p.seq[length(p.seq)])
  })

  reset.rownames(
    do.call('rbind.data.frame', new.p.start.end)
  )
}
