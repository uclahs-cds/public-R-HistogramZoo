
identify.uniform.segments = function(
  seg.gr,
  bin.counts,
  trim.peak.threshold,
  short.peak.threshold
){

  # Peak width and threshold
  peak.width = GenomicRanges::width(seg.gr)
  unif.threshold = floor((1-trim.peak.threshold)*peak.width)
  if(short.peak.threshold < 1){
    short.peak.threshold = short.peak.threshold*peak.width
  } else {
    short.peak.threshold = rep(short.peak.threshold, length(peak.width))
  }

  # Initializing Variables & GRanges
  chr.gene = GenomicRanges::seqnames(seg.gr)[1]
  strand.gene = GenomicRanges::strand(seg.gr)[1]
  unif.seg.gr = GenomicRanges::GRanges()

  # Peak trimming
  for(i in 1:length(seg.gr)){

    # Initializing
    seg.start = GenomicRanges::start(seg.gr)[i]
    seg.end = GenomicRanges::end(seg.gr)[i]
    coverage.segment = bin.counts$Coverage[bin.counts$start >= seg.start & bin.counts$end <= seg.end]

    # Looking for uniform stretches below the threshold
    seg.pattern = diff(coverage.segment)
    seg.breakpoints = c(0, which(seg.pattern != 0), nrow(coverage.segment))
    zero.stretch = diff(seg.breakpoints)
    start.seg.bp = seg.breakpoints[which(zero.stretch >= unif.threshold[i])]+seg.start
    end.seg.bp = seg.breakpoints[which(zero.stretch >= unif.threshold[i])+1]+seg.start-1

    # Segmenting a peak
    if(length(start.seg.bp) > 0){
      unif.gr = GenomicRanges::GRanges(seqnames = chr.gene,
                                       IRanges::IRanges(start = start.seg.bp, end = end.seg.bp),
                                       strand = strand.gene)
      nonunif.gr = GenomicRanges::setdiff(seg.gr[i], unif.gr)
      peak.seg.gr = c(unif.gr, nonunif.gr)
      peak.seg.gr = peak.seg.gr[GenomicRanges::width(peak.seg.gr) >= short.peak.threshold[i]]
    } else { # No uniform segments are found
      peak.seg.gr = seg.gr[i]
    }
    peak.seg.gr$i = i

    # Add to master GRanges object
    unif.seg.gr = c(unif.seg.gr, peak.seg.gr)
    unif.seg.gr = sort(unif.seg.gr)
  }

  return(unif.seg.gr)
}
