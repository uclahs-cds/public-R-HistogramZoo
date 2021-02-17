.generate.null.result = function(PARAMETERS){

  # Required Column Names
  bed.cols = c("chr", "start", "end", "name", "score", "strand",
               "thickStart", "thickEnd", "itemRgb",
               "blockCount", "blockSizes","blockStarts")

  # Number of Columns
  output.ncol = length(PARAMETERS$ALL.SAMPLES) + 12

  # A table with no values
  output.table = data.frame(matrix(ncol = output.ncol, nrow = 0))
  names(output.table) = c(bed.cols, PARAMETERS$ALL.SAMPLES)

  return(output.table)
}
