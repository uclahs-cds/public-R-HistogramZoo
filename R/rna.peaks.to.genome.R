

.rna.peaks.to.genome = function(merged.peaks.rna, GENEINFO){
  
  # Return an empty GRanges if no peaks survive
  if(length(merged.peaks.rna) == 0){
    return(GenomicRanges::GRanges())
  }
  
  # Transferring to Genomic Coordinates
  merged.peaks.genome = merged.peaks.rna
  GenomicRanges::end(merged.peaks.genome) = GENEINFO$RNA2DNA[GenomicRanges::end(merged.peaks.genome)]
  GenomicRanges::start(merged.peaks.genome) = GENEINFO$RNA2DNA[GenomicRanges::start(merged.peaks.genome)]
  
  # Creating a GRanges object from the Annotation
  anno_gr = GenomicRanges::makeGRangesFromDataFrame(GENEINFO$anno)
  
  # Filtering Out Introns
  merged.peaks.filtered.genome = GenomicRanges::GRanges()
  for(i in 1:length(merged.peaks.genome)){
    tmp.gr = GenomicRanges::intersect(merged.peaks.genome[i], anno_gr)
    if(length(tmp.gr) > 0){S4Vectors::mcols(tmp.gr) = S4Vectors::mcols(merged.peaks.genome[i])}
    merged.peaks.filtered.genome = c(merged.peaks.filtered.genome, tmp.gr)
  }
  
  return(merged.peaks.filtered.genome)
}