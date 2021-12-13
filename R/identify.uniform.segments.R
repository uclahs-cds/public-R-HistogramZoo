
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

# meaningful.interval(
#   h = x/sum(x),
#   p = generate.unif(x),
#   a = todo$Var1[i],
#   b = todo$Var2[i],
#   N = sum(x),
#   L = length(x)
# )

#' Finds the largest uniform segment that is longer than threshold
find.uniform.segment = function(x, metric = c("jaccard", "intersection", "ks", "mse", "chisq"), threshold = 0.5, step.size = 1, max.sd.size = 1) {
  metric = match.arg(metric)
  num.bins = length(x)
  min.seg.size = ceiling(num.bins * threshold)
  metric.func = get(paste('histogram', metric, sep = "."))

  p.unif = generate.unif(x)
  res = lapply(seq(from = 1, to = num.bins - min.seg.size, by = step.size), function(a) {
    lapply(seq(from = min.seg.size + a, to = num.bins, by = step.size), function(b) {
      x.sub = x[a:b]
      p.unif.sub = generate.unif(x.sub)
      h.sub = x.sub / sum(x.sub)

      m = metric.func(h.sub, p.unif.sub)

      list(a = a, b = b, metric = m)
    })
  })

  res.df = do.call(rbind.data.frame, unlist(res, recursive = F))
  colnames(res.df) = c("a", "b", "metric")
  res.df$length = res.df$b - res.df$a

  # Select the longest interval that is within 1 sd of the maximum
  min.metric = min(res.df$metric, na.rm = T)
  sd.metric = sd(res.df$metric, na.rm = T)
  sd.metric = ifelse(is.na(sd.metric), 0, sd.metric) # If there's only 1 case
  # The range in which we are looking for the minimum
  res.sd.range = res.df[res.df$metric <= min.metric + sd.metric * max.sd.size, ]
  min.interval.index = which.min(res.sd.range$length)

  as.list(res.sd.range[min.interval.index,])
}

# Uniform Generation
generate.unif = function(x){
  rep(1/length(x), length(x))
}
