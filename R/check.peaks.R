.check.peaks = function(PEAKS){
  peak.names = c("chr", "start", "end", "name", "score", "strand", "blockCount", "blockSizes", "blockStarts", "sample")

  if(!all(peak.names %in% colnames(PEAKS))){stop("Missing Column(s) in PEAKS")}
  if(!is.character(PEAKS$chr)){stop("The chr column in PEAKS neads to be character")}
  if(!is.integer(PEAKS$start) | !is.integer(PEAKS$end)){stop("The start and end columns in PEAKS have to be integer")}
  if(!is.character(PEAKS$name)){stop("The name column in PEAKS needs to be character")}
  if(!all(PEAKS$strand %in% c("+", "-", "*"))){stop("The strand column must be one of +, -, *")}
  if(!is.numeric(PEAKS$score)){stop("The score column in PEAKS needs to be numeric")}
  if(!is.character(PEAKS$sample)){stop("The sample column in PEAKS must be character")}

  # Can check that blockCount, blockSizes and blockStarts actually fulfill correct formatting
}
