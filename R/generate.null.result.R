.generate.null.result = function(all.samples){

  # Required Column Names
  bed.cols = c("chr", "start", "end", "name", "score", "strand",
               "thickStart", "thickEnd", "itemRgb",
               "blockCount", "blockSizes","blockStarts")

  # Number of Columns
  output.ncol = length(all.samples) + 12

  # A table with no values
  output.table = data.frame(matrix(ncol = output.ncol, nrow = 0))
  names(output.table) = c(bed.cols, all.samples)

  return(output.table)
}
