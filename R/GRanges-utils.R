#' Converts a GRanges object from base 0 to base 1
#'
#' @param gr a GRanges object in base 0
#' @return a GRanges object in base 1
base0_to_base1 <- function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)+1
  return(gr)
}

#' Generates region_ids for GRanges objects
#'
#' @param gr a GRanges object
#' @return character identifiers for ranges in a GRanges object
generate_identifiers <- function(gr){
  return(
    paste0(
      GenomicRanges::seqnames(gr), ":",
      GenomicRanges::start(gr), "-",
      GenomicRanges::end(gr), ":",
      GenomicRanges::strand(gr)
    )
  )
}

#' Generates region_ids for GRangesList objects
#'
#' @param gr a GRangesList object
#' @return character identifiers for each GRanges object in a GRangesList
generate_grangeslist_identifiers <- function(gr_list){
  gr_ranges <- range(gr_list)
  return(
    paste0(
      GenomicRanges::seqnames(gr_ranges), ":",
      GenomicRanges::start(gr_ranges), "-",
      GenomicRanges::end(gr_ranges), ":",
      GenomicRanges::strand(gr_ranges)
    )
  )
}
