.merge.p = function(PEAKSGR, MERGED.PEAKS, ANNOTATION, PARAMETERS, ID.COLS){

  # Numeric Stability
  PEAKSGR$score = ifelse(PEAKSGR$score == 0, 2e-16, PEAKSGR$score)

  # Tag
  MERGED.PEAKS$tag = apply(GenomicRanges::mcols(MERGED.PEAKS)[,ID.COLS], 1, function(x) paste(x, collapse = ":"))

  # Finding Overlaps
  overlaps = GenomicRanges::findOverlaps(MERGED.PEAKS, PEAKSGR)
  overlapping_peaks = PEAKSGR[S4Vectors::subjectHits(overlaps)]
  overlapping_peaks$merged_peak = MERGED.PEAKS$tag[S4Vectors::queryHits(overlaps)]
  overlapping_peaks$tag = paste0(overlapping_peaks$sample, ":", overlapping_peaks$merged_peak)

  overlapping_peaks = split(overlapping_peaks, overlapping_peaks$tag)
  results = do.call(rbind, lapply(1:length(overlapping_peaks), function(i){
    tmp = overlapping_peaks[[i]]
    unique_p = unique(GenomicRanges::mcols(tmp)[c("peak", "score")])
    p = ifelse(nrow(unique_p) > 1, metap::sumlog(unique_p$score)$p, unique_p$score)
    data.frame(
      "peak" = tmp$merged_peak[1],
      "sample" = tmp$sample[1],
      "p" = p,
      stringsAsFactors = F
    )
  }))

  # Format this into something usable
  pmat = reshape2::dcast(results, peak ~ sample, value.var = "p")
  pmat[is.na(pmat)] <- 1

  # Adding Extra Samples
  missing.samples = setdiff(PARAMETERS$ALL.SAMPLES, colnames(pmat))
  if(length(missing.samples) > 0){
    add.table = data.frame(matrix(1, nrow = nrow(pmat), ncol = length(missing.samples), dimnames = list(NULL, missing.samples)))
    pmat = cbind(pmat, add.table)
  }
  pmat = pmat[,c("peak", PARAMETERS$ALL.SAMPLES), drop = FALSE]

  return(pmat)

}
