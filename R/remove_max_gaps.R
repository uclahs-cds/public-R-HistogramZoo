#' Creates a new set of segments from a partition of points
#'
#' @param start_end_points a data.frame with start and end named columns representing intervals
#' @param max_gaps a data.frame with start and end named columns representing gaps
#' @param remove_short_segment the minimum segment length to include in the final segments
#'
#' @return a data.frame of the intervals with max gaps removed
#' \describe{
#'     \item{start}{start indices of new intervals}
#'     \item{end}{end indices of new intervals}
#' }
#' @export
remove_max_gaps <- function(start_end_points, max_gaps, remove_short_segment = 0) {

  # Error checking
  if(!is.data.frame(start_end_points) | !all(c("start", "end") %in% colnames(start_end_points))){
    stop("start_end_points should be a data.frame with columns start and end")
  }
  if(!is_equal_integer(start_end_points[,"start"]) | !is_equal_integer(start_end_points[,"end"]) | !all(start_end_points[,"start"] <= start_end_points[,"end"])){
    stop("start_end_points start and end columns must represent functional intervals, i.e end >= start and end and start are numeric integers")
  }
  if(nrow(max_gaps) == 0) {
    # Return the original points as start/end data frame
    return(start_end_points)
  }
  if(!is.data.frame(max_gaps) | !all(c("start", "end") %in% colnames(max_gaps))){
    stop("max_gaps should be a data.frame with columns start and end")
  }
  if(!is_equal_integer(max_gaps[,"start"]) | !is_equal_integer(max_gaps[,"end"]) | !all(max_gaps[,"start"] <= max_gaps[,"end"])){
    stop("max_gaps start and end columns must represent functional intervals, i.e end >= start and end and start are numeric integers")
  }
  if(!is_equal_integer(remove_short_segment) | length(remove_short_segment) != 1){
    stop("remove_short_segment must be a numeric integer of length 1")
  }

  # For each segment, create a sequence of consecutive integers
  seg_seq_list <- mapply(seq.int, from = start_end_points$start, to = start_end_points$end, by = 1, SIMPLIFY = FALSE)

  max_gaps_seq <- unlist(lapply(1:nrow(max_gaps), function(i) {
    # TODO: Should this include the endpoints or no?
    seq(max_gaps[i, 1], max_gaps[i, 2], by = 1)
  }))

  seg_seq_no_maxgap <- lapply(seg_seq_list, function(p_seq) {
    sort(unique(setdiff(p_seq, max_gaps_seq)))
  })

  # https://stackoverflow.com/a/24837419
  # Take the segments and regroup into consecutive integers
  new_p_seq <- lapply(seg_seq_no_maxgap, function(p_no_maxgap) {
    consec_p <- split(p_no_maxgap, cumsum(c(1, diff(p_no_maxgap) != 1)))
    # Remove short segments
    # TODO: Should we do this later? Do we need to join any segments first?
    consec_p[unlist(lapply(consec_p, length) > remove_short_segment)]
  })

  # Remove one level of the lists for a flat structure
  new_p_seq_consec <- unname(unlist(new_p_seq, recursive = FALSE))

  # Extract the start and end points from the sequences
  new_p_start_end <- lapply(new_p_seq_consec, function(p_seq) {
    list('start' = p_seq[1], 'end' = p_seq[length(p_seq)])
  })

  res <- do.call('rbind.data.frame', new_p_start_end)
  rownames(res) <- NULL

  return(res)
}
