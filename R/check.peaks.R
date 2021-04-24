.check.peaks = function(peaks){
  peak.names = c("chr", "start", "end", "name", "score", "strand", "blockCount", "blockSizes", "blockStarts", "sample")

  if(!all(peak.names %in% colnames(peaks))){stop("Missing Column(s) in peaks")}
  if(!is.character(peaks$chr)){stop("The chr column in peaks neads to be character")}
  if(!is.integer(peaks$start) | !is.integer(peaks$end)){stop("The start and end columns in peaks have to be integer")}
  if(!is.character(peaks$name)){stop("The name column in peaks needs to be character")}
  if(!all(peaks$strand %in% c("+", "-", "*"))){stop("The strand column must be one of +, -, *")}
  if(!is.numeric(peaks$score)){stop("The score column in peaks needs to be numeric")}
  if(!is.character(peaks$sample)){stop("The sample column in peaks must be character")}

  # Can check that blockCount, blockSizes and blockStarts actually fulfill correct formatting
}
