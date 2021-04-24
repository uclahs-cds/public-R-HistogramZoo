.bed6tobed12 = function(
  merged.peaks,
  id.cols
){
  if(length(id.cols) == 1){
    merged.peaks$tag = GenomicRanges::mcols(merged.peaks)[,id.cols]
  } else{
    merged.peaks$tag = apply(GenomicRanges::mcols(merged.peaks)[,id.cols], 1, function(x) paste(x, collapse = ":"))
  }

  bed12 = do.call(rbind, lapply(unique(merged.peaks$tag), function(itag) {
    tmp.peak = merged.peaks[merged.peaks$tag == itag]
    tmp.peak = tmp.peak[order(GenomicRanges::start(tmp.peak))]
    data.frame(
      "chr" = as.character(GenomicRanges::seqnames(tmp.peak))[1],
      "start" = min(GenomicRanges::start(tmp.peak)),
      "end" = max(GenomicRanges::end(tmp.peak)),
      "name" = GenomicRanges::mcols(tmp.peak)$name[1],
      "score" = 0, # This might be changed to a p-value if there is one
      "strand" = as.character(GenomicRanges::strand(tmp.peak))[1],
      "thickStart" = min(GenomicRanges::start(tmp.peak)),
      "thickEnd" = max(GenomicRanges::end(tmp.peak)),
      "itemRgb" = 0,
      "blockCount" = length(tmp.peak),
      "blockSizes" = paste0(paste(GenomicRanges::end(tmp.peak) - GenomicRanges::start(tmp.peak), collapse = ","), ","),
      "blockStarts" = paste0(paste(GenomicRanges::start(tmp.peak) - min(GenomicRanges::start(tmp.peak)), collapse = ","), ","),
      "peak" = itag,
      stringsAsFactors = F
      )
  }))
  return(bed12)
}
