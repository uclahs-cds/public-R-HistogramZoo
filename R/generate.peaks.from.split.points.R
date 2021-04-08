
generate.peaks.from.split.points = function(
  p,
  GENEPEAKSGR,
  GENEINFO,
  m = 100
) {
  
  # Generate Regions With Peak Coverage
  REDUCED.GENE.PEAKS.GR = GenomicRanges::reduce(GENEPEAKSGR)
  
  # Find Points Within Peaks
  POINT.GR = GenomicRanges::GRanges(seqnames = GENEINFO$chr, IRanges::IRanges(p), strand = GENEINFO$strand)
  POINT.GR = POINT.GR[S4Vectors::queryHits(GenomicRanges::findOverlaps(POINT.GR, REDUCED.GENE.PEAKS.GR))]
  
  # Remove Split Points Close To Peak Boundaries
  ovl = GenomicRanges::findOverlaps(POINT.GR, REDUCED.GENE.PEAKS.GR)
  dist.to.start = (GenomicRanges::start(POINT.GR)[S4Vectors::queryHits(ovl)] - 
                     GenomicRanges::start(REDUCED.GENE.PEAKS.GR)[S4Vectors::subjectHits(ovl)])
  dist.to.end = (GenomicRanges::end(REDUCED.GENE.PEAKS.GR)[S4Vectors::subjectHits(ovl)] - 
                   GenomicRanges::start(POINT.GR)[S4Vectors::queryHits(ovl)])
  POINT.GR = POINT.GR[dist.to.start > m & dist.to.end > m,]
  
  # Return Peaks Object If All Split Points Fall Outside of Peak Regions 
  if(length(POINT.GR) == 0){
    return(REDUCED.GENE.PEAKS.GR)
  }
  
  # Split Segments
  p.midpoints = GenomicRanges::start(POINT.GR)
  p.all = c(p.midpoints, GenomicRanges::start(REDUCED.GENE.PEAKS.GR), GenomicRanges::end(REDUCED.GENE.PEAKS.GR))
  p.all = sort(p.all)
  SEG.GR = GenomicRanges::GRanges(
    seqnames = GENEINFO$chr,
    IRanges::IRanges(
      start = p.all[1:(length(p.all)-1)], 
      end = p.all[2:length(p.all)]), 
    strand = GENEINFO$strand)
  
  # Shifting Segments by 1 bp if it's a midpoint
  GenomicRanges::start(SEG.GR) = ifelse(GenomicRanges::start(SEG.GR) %in% p.midpoints, GenomicRanges::start(SEG.GR) + 1, GenomicRanges::start(SEG.GR))

  return(SEG.GR)
  
}