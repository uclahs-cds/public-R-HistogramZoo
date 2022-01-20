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
