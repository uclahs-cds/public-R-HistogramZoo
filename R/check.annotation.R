.check.annotation <- function(annotation, peaksgr, gene){

  # Checking that all columns exist
  annot.cols = c("chr", "feature", "start", "stop", "strand", "gene", "transcript")
  if(!all(annot.cols %in% names(annotation))){return(FALSE)}

  # Subsetting the annotation into the ranges of the gene
  annotation = annotation[annotation$gene == gene,]
  annot.gr = GenomicRanges::makeGRangesFromDataFrame(annotation, keep.extra.columns = T)

  # Converting to base 1
  GenomicRanges::start(peaksgr) = GenomicRanges::start(peaksgr)+1

  # Checking to see that all peaks fall within the annotation
  check.annot = GenomicRanges::setdiff(peaksgr, annot.gr)

  if(length(check.annot) == 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
