remove.max.gaps.agnostic = function(p, max.gaps, remove.short.segment = 1) {
  p.seq = unlist(lapply(2:length(p), function(i) {
    seq(p[i - 1], p[i], by = 1)
    }))

  max.gaps.seq = unlist(lapply(1:nrow(max.gaps), function(i) {
    seq(max.gaps[i, 1], max.gaps[i, 2], by = 1)
  }))

  p.no.maxgap = sort(unique(setdiff(p.seq, max.gaps.seq)))

  # https://stackoverflow.com/a/24837419
  # Take the segments and regroup into consecutive integers
  new.p.seq <- split(p.no.maxgap, cumsum(c(1, diff(p.no.maxgap) != 1)))
  new.p.seq <- new.p.seq[unlist(lapply(new.p.seq, length) > remove.short.segment)]

  # Only take the first and last consecutive numbers
  new.p <- lapply(new.p.seq, function(y) {
    list(start = y[1], end = y[length(y)])
  })

  new.p
}

remove.max.gaps.iranges = function(p, max.gaps, remove.short.segment = 1){
  # Create IRanges object with the segments identified by p
  l.p = length(p) - 1
  p.ir = IRanges::IRanges(start = p[1:l.p], end = c(p[2:l.p] - 1, p[l.p+1]))
  p.ir = S4Vectors::split(p.ir, seq_along(p.ir))
  # Create IRanges object with max gap
  maxgap.ir = if(nrow(max.gaps) > 0) IRanges::IRanges(start = max.gaps[,1], end = max.gaps[,2]) else IRanges::IRanges()
  # Subtract maxgaps
  segs.ir = IRanges::setdiff(p.ir, maxgap.ir)
  segs.ir = unlist(segs.ir)
  # Remove short segments
  segs.ir = segs.ir[IRanges::width(segs.ir) > remove.short.segment]
  # Create a list of ranges
  lapply(seq_along(segs.ir), function(i) list(start = IRanges::start(segs.ir[i]), end = IRanges::end(segs.ir[i])))
}
