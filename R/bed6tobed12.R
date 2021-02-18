.bed6tobed12 = function(
  MERGED.PEAKS,
  ID.COLS
){
  if(length(ID.COLS) == 1){
    MERGED.PEAKS$tag = GenomicRanges::mcols(MERGED.PEAKS)[,ID.COLS]
  } else{
    MERGED.PEAKS$tag = apply(GenomicRanges::mcols(MERGED.PEAKS)[,ID.COLS], 1, function(x) paste(x, collapse = ":"))
  }

  BED12 = do.call(rbind, lapply(unique(MERGED.PEAKS$tag), function(itag) {
    TMP.PEAK = MERGED.PEAKS[MERGED.PEAKS$tag == itag]
    TMP.PEAK = TMP.PEAK[order(GenomicRanges::start(TMP.PEAK))]
    bed12 = data.frame(
      "chr" = as.character(GenomicRanges::seqnames(TMP.PEAK))[1],
      "start" = min(GenomicRanges::start(TMP.PEAK)),
      "end" = max(GenomicRanges::end(TMP.PEAK)),
      "name" = GenomicRanges::mcols(TMP.PEAK)$name[1],
      "score" = 0, # This might be changed to a p-value if there is one
      "strand" = as.character(GenomicRanges::strand(TMP.PEAK))[1],
      "thickStart" = min(GenomicRanges::start(TMP.PEAK)),
      "thickEnd" = max(GenomicRanges::end(TMP.PEAK)),
      "itemRgb" = 0,
      "blockCount" = length(TMP.PEAK),
      "blockSizes" = paste0(paste(GenomicRanges::end(TMP.PEAK) - GenomicRanges::start(TMP.PEAK), collapse = ","), ","),
      "blockStarts" = paste0(paste(GenomicRanges::start(TMP.PEAK) - min(GenomicRanges::start(TMP.PEAK)), collapse = ","), ","),
      "peak" = itag,
      stringsAsFactors = F
      )

  }))
  return(BED12)
}
