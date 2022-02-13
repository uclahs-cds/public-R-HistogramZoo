base0.to.base1 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)+1
  return(gr)
}

base1.to.base0 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)-1
  return(gr)
}

generate.identifiers = function(gr){
  paste0(GenomicRanges::seqnames(gr), ":", GenomicRanges::start(gr), "-", GenomicRanges::end(gr))
}
