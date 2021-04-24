.retrieve.peaks.as.granges = function(peaks, gene, return.df = F){

  peaks = peaks[peaks$name == gene,]

  # Creating a bed6 File
  bed6 = data.frame(stringsAsFactors = F)
  for(i in 1:nrow(peaks)){
    tmp = data.frame(
      "chr" = peaks$chr[i],
      "start" = peaks$start[i] + strtoi(strsplit(peaks$blockStarts[i], split = ",")[[1]]),
      "end" = peaks$start[i] + strtoi(strsplit(peaks$blockStarts[i], split = ",")[[1]]) + strtoi(strsplit(peaks$blockSizes[i], split = ",")[[1]]),
      "name" = peaks$name[i],
      "score" = peaks$score[i],
      "strand" = peaks$strand[i],
      "peak" = i,
      "sample" = peaks$sample[i],
      stringsAsFactors = F
    )
    bed6 = rbind(bed6, tmp)
  }
  if(return.df){
    return(bed6)
  } else {
    return(GenomicRanges::makeGRangesFromDataFrame(bed6, keep.extra.columns = T))
  }
}
