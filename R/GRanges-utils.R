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
