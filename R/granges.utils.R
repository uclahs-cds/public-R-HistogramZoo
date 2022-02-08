
split.strand = function(gr){
  res = S4Vectors::split(gr, f = GenomicRanges::strand(gr))
  res[sapply(res, length) > 0]
}

split.id = function(gr, reduce = T){
  if(is(gr, "GRanges")){
    res = S4Vectors::split(gr, f = gr$id)
    if(reduce) {res = GenomicRanges::reduce(res)}
  } else if (is(gr, "CompressedGRangesList")) {
    res = lapply(gr, function(x) S4Vectors::split(x, x$id))
    if(reduce) {
      res = lapply(res, GenomicRanges::reduce)
    }
  }
  res
}

ApplyGRangesFunc = function(gr, func, ...){
  if(is(gr, "GRanges")){
    func(gr, ...)
  } else if( is(gr, "CompressedGRangesList") ){
    lapply(gr, function(i) func(i, ...))
  }
}

# May or may not need this function
filter.regions = function(regions.gr, by.gr){
  if(is(regions.gr, "CompressedGRangesList")){
    select.regions = S4Vectors::queryHits(GenomicRanges::findOverlaps(regions.gr, by.gr))
    res = regions.gr[sort(unique(select.regions))]
  } else if (is(regions.gr, "list") & all(sapply(regions.gr, function(x) is(x, "CompressedGRangesList")))) {
    res = lapply(regions.gr, function(reg.gr){
      select.regions = S4Vectors::queryHits(GenomicRanges::findOverlaps(regions.gr, by.gr))
      regions.gr[sort(unique(select.regions))]
    })
  }
  res
}

base0.to.base1 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)+1
  return(gr)
}

base1.to.base0 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)-1
  return(gr)
}
