.bed6tobed12 = function(
  merged.peaks,
  name.id,
  peak.id,
  score.id,
  distribution.id,
  meta.id
){

  merged.peaks$tag = paste0(GenomicRanges::mcols(merged.peaks)[,name.id], ":", GenomicRanges::mcols(merged.peaks)[,peak.id])

  bed12 = do.call(rbind, lapply(unique(merged.peaks$tag), function(itag) {
    tmp.peak = merged.peaks[merged.peaks$tag == itag]
    tmp.peak = tmp.peak[order(GenomicRanges::start(tmp.peak))]
    data.frame(
      "chr" = as.character(GenomicRanges::seqnames(tmp.peak))[1],
      "start" = min(GenomicRanges::start(tmp.peak)),
      "end" = max(GenomicRanges::end(tmp.peak)),
      "name" = GenomicRanges::mcols(tmp.peak)$name[1],
      "score" = GenomicRanges::mcols(tmp.peak)[1, score.id],
      "strand" = as.character(GenomicRanges::strand(tmp.peak))[1],
      "thickStart" = min(GenomicRanges::start(tmp.peak)),
      "thickEnd" = max(GenomicRanges::end(tmp.peak)),
      "itemRgb" = 0,
      "blockCount" = length(tmp.peak),
      "blockSizes" = paste0(paste(GenomicRanges::end(tmp.peak) - GenomicRanges::start(tmp.peak), collapse = ","), ","),
      "blockStarts" = paste0(paste(GenomicRanges::start(tmp.peak) - min(GenomicRanges::start(tmp.peak)), collapse = ","), ","),
      "distribution" = GenomicRanges::mcols(tmp.peak)[1, distribution.id],
      "distribution.features" = GenomicRanges::mcols(tmp.peak)[1, meta.id],
      "peak" = itag,
      stringsAsFactors = F
      )
  }))
  return(bed12)
}
