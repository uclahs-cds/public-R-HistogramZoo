#' Converts a GRanges object from base 0 to base 1
#' @param gr A GRanges object
base0_to_base1 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)+1
  return(gr)
}

#' Converts a GRanges object from base 1 to base 0
#' @param gr A GRanges object
base1_to_base0 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)-1
  return(gr)
}

#' Generates region_ids for GRanges objects
#' @param gr A GRanges object
generate_identifiers = function(gr){
  return(
    paste0(
      GenomicRanges::seqnames(gr), ":", 
      GenomicRanges::start(gr), "-", 
      GenomicRanges::end(gr), ":", 
      GenomicRanges::strand(gr)
    ) 
  )
}
