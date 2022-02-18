
# Uniform Generation
generate.unif = function(x){
  rep(1/length(x), length(x))
}

calc.prob.diff = function(h, p, a, b){
  interval = a:b
  # Round to prevent floating point issues
  hab = round(sum(h[interval]), digits = 14)
  pab = round(sum(p[interval]), digits = 14)
  hab > pab
}

# Meaningful interval
meaningful.interval = function(h, p, a, b, N, L){
  relative.entropy = rel.entropy(h, p, a, b)
  prob.diff = calc.prob.diff(h, p, a, b)
  c(mint = relative.entropy >= (1/N)*log(L*(L+1)/2) && prob.diff, entropy = relative.entropy)
}

# Meaningful gap
meaningful.gap = function(h, p, a, b, N, L){
  relative.entropy = rel.entropy(h, p, a, b)
  if(is.na(relative.entropy)) relative.entropy = Inf
  prob.diff = calc.prob.diff(h, p, a, b)
  mgap = (relative.entropy >= (1/N)*log(L*(L+1)/2) && !prob.diff) || (all(h == 0))
  c(mgap = mgap, entropy = relative.entropy)
}

maximal.meaningful = function(x) {
  curr.df = x
  max.intervals = data.frame()
  while(nrow(curr.df) > 0) {
    maximum.entropy.index = which.max(curr.df$entropy)
    maximum.entropy = curr.df[maximum.entropy.index, ]
    maximum.entropy.seq = seq(maximum.entropy$start, maximum.entropy$end)
    max.intervals = rbind(max.intervals, maximum.entropy)

    # Find all of the segments that overlap.
    # These will all be less than the maximum
    overlap.max = mapply(function(from, to) {
      s = seq(from, to)
      length(intersect(maximum.entropy.seq, s)) > 0
    }, from = curr.df$start, to = curr.df$end)

    curr.df = curr.df[!overlap.max, ]
  }
  # Preserve the old index
  max.intervals$index = rownames(max.intervals)
  rownames(max.intervals) <- NULL
  max.intervals
}


#' Find all meaningful gaps
#'
#' @param x histogram (vector of counts)
#' @param change.points Change points
#'
#' @return TODO
#' @export
find.all.meaningful.gap = function(x, change.points) {
  todo = expand.grid(start = change.points, end = change.points)
  todo = todo[todo$end > todo$start,]

  mgap = do.call(rbind, lapply(1:nrow(todo), function(i) {
    meaningful.gap(
      h = x/sum(x),
      p = generate.unif(x),
      a = todo$start[i],
      b = todo$end[i],
      N = sum(x),
      L = length(x)
    )
  }))
  df = cbind(todo, mgap)
  # df = df[order(df$end, df$start),]
  df = df[order(df$entropy), ]
  df$mgap = as.numeric(df$mgap)
  # df$scaled_entropy = (df$entropy - min(df$entropy, na.rm = T)) / (max(df$entropy, na.rm = T) - min(df$entropy, na.rm = T))

  seg.gap.data = df[df$mgap > 0 & !is.na(df$mgap), ]
  # seg.gap.data$index = 1:nrow(seg.gap.data)

  maximal.meaningful(seg.gap.data)
}

#' Finds the meaningful gaps between the points in s
#'
#' @param x The histogram data
#' @param seg.points the segment points
#' @param change.points the change points
#' @param min.gap The minimum gap to be considered a meaningful gap
#'
#' @return TODO
#' @export
meaningful.gaps.local = function(x, seg.points, change.points, min.gap = 2) {

  max.gaps.list <- lapply(seq(2, length(seg.points)), function(i) {
    x.sub = x[seg.points[i-1]:seg.points[i]]
    chg.pts = change.points[change.points >= seg.points[i-1] & change.points <= seg.points[i]] - seg.points[i-1] + 1

    max.gaps = find.all.meaningful.gap(x.sub, chg.pts)

    if(nrow(max.gaps) > 0) {
      max.gaps[, c('start','end')] <- max.gaps[, c('start','end')] + seg.points[i-1] - 1
      max.gaps$seg.start = seg.points[i-1]
      max.gaps$seg.end = seg.points[i]
      max.gaps
    }
  })

  max.gaps <- do.call(rbind.data.frame, max.gaps.list)
  # Remove gaps that are smaller than min.gap
  max.gaps[max.gaps$end - max.gaps$start >= min.gap, ]
}

find.new.segments = function(gaps.df) {
  new.segments = c()
  for(i in rownames(gaps.df)) {
    r = gaps.df[i, ]
    # Have a large gap
    if(r$end - r$start > 1) {
      left.diff = abs(r$seg.start - r$start)
      right.diff = abs(r$seg.end-r$end)
      if(left.diff < right.diff && right.diff > 2) {
        new.segments = c(new.segments, r$end - 1)
      } else if(left.diff > right.diff && left.diff > 2){
        new.segments = c(new.segments, r$start + 1)
      }
    }
  }
  new.segments
}
