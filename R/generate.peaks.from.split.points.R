
generate.peaks.from.split.points = function(
  p,
  genepeaksgr,
  geneinfo,
  m = 100
) {

  # Generate Regions With Peak Coverage
  reduced.gene.peaks.gr = GenomicRanges::reduce(genepeaksgr)

  # If there are no fitted points
  if(length(p) == 0){
    return(reduced.gene.peaks.gr)
  }

  # Find Points Within Peaks
  point.gr = GenomicRanges::GRanges(seqnames = geneinfo$chr, IRanges::IRanges(p), strand = geneinfo$strand)
  point.gr = point.gr[S4Vectors::queryHits(GenomicRanges::findOverlaps(point.gr, reduced.gene.peaks.gr))]

  # Remove Split Points Close To Peak Boundaries
  ovl = GenomicRanges::findOverlaps(point.gr, reduced.gene.peaks.gr)
  dist.to.start = (GenomicRanges::start(point.gr)[S4Vectors::queryHits(ovl)] -
                     GenomicRanges::start(reduced.gene.peaks.gr)[S4Vectors::subjectHits(ovl)])
  dist.to.end = (GenomicRanges::end(reduced.gene.peaks.gr)[S4Vectors::subjectHits(ovl)] -
                   GenomicRanges::start(point.gr)[S4Vectors::queryHits(ovl)])
  point.gr = point.gr[dist.to.start > m & dist.to.end > m,]

  # Return Peaks Object If All Split Points Fall Outside of Peak Regions
  if(length(point.gr) == 0){
    return(reduced.gene.peaks.gr)
  }

  ##################
  # Split Segments #
  ##################

  # Setting point boundaries
  p.start = GenomicRanges::start(reduced.gene.peaks.gr)
  p.end = GenomicRanges::end(reduced.gene.peaks.gr)
  p.midpoints = GenomicRanges::start(point.gr)
  p.all = c(p.midpoints, p.start, p.end)
  p.all = sort(p.all)

  # Initializing segments
  seg.gr = GenomicRanges::GRanges(
    seqnames = geneinfo$chr,
    IRanges::IRanges(
      start = p.all[1:(length(p.all)-1)],
      end = p.all[2:length(p.all)]),
    strand = geneinfo$strand)

  # Shifting Segments by 1 bp if it's a segment point or an endpoint of a peak
  GenomicRanges::start(seg.gr) = ifelse(GenomicRanges::start(seg.gr) %in% c(p.midpoints, p.end), GenomicRanges::start(seg.gr) + 1, GenomicRanges::start(seg.gr))
  GenomicRanges::end(seg.gr) = ifelse(GenomicRanges::end(seg.gr) %in% c(p.start), GenomicRanges::end(seg.gr) - 1, GenomicRanges::end(seg.gr))

  # Remove the Non-peak segments
  seg.gr = seg.gr[unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(seg.gr, reduced.gene.peaks.gr)))]

  return(seg.gr)

}
