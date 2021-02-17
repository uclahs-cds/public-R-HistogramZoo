
.check.annotation <- function(ANNOTATION, PEAKSGR, GENE){

  # Checking that all columns exist
  annot.cols = c("chr", "feature", "start", "stop", "strand", "gene", "transcript")
  if(!all(annot.cols %in% names(ANNOTATION))){
    return(FALSE)
  }
  
  # Subsetting the annotation into the ranges of the gene
  ANNOTATION = ANNOTATION[ANNOTATION$gene == GENE,]
  annot.gr = GenomicRanges::makeGRangesFromDataFrame(ANNOTATION, keep.extra.columns = T)
  
  # Converting to base 1
  GenomicRanges::start(PEAKSGR) = GenomicRanges::start(PEAKSGR)+1
  
  # Checking to see that all peaks fall within the annotation
  check.annot = GenomicRanges::setdiff(PEAKSGR, annot.gr)
  
  if(length(check.annot) == 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
