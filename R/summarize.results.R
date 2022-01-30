#' Formats results
#'
#' @param coverage.model.obj
#' @param output.format
#'
#' @return
#' @export
#'
#' @examples
#' TODO: Create a reasonable format for basic histograms
summarize.results = function(
  coverage.model.obj,
  output.format = c("stats.only", "BED6", "BED12")
){

  output.format = match.arg(output.format)

  # Coverage.model.obj
  gene.model = coverage.model.obj$gene.model
  histogram.ids = names(gene.model)
  histogram.bin.size = coverage.model.obj$histogram.bin.size
  results = lapply(coverage.model.obj$results, `[[`, "models")

  formatted.results = GenomicRanges::GRanges()
  for(i in histogram.ids){
    x = gene.model[i]
    x = unlist(x)
    bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
    bins = unlist(bins)
    res = results[[i]]
    for(j in seq_along(res)){
      id =  paste0(i, ":", j)
      final.mod = res[[j]]$majority.vote
      f.res = GenomicRanges::reduce(bins[final.mod$seg.start:final.mod$seg.end])
      f.res = base1.to.base0(f.res)
      f.res = if(output.format == "BED12") bed6tobed12(f.res)
      S4Vectors::mcols(f.res)$name = id
      S4Vectors::mcols(f.res)[,c( "value", "metric", "dist")] = final.mod[c( "value", "metric", "dist")]
      S4Vectors::mcols(f.res)$params = dput.str(final.mod$par)
      formatted.results[id] = f.res
    }
  }
  # Reordering the columns
  S4Vectors::mcols(formatted.results) = S4Vectors::mcols(formatted.results)[,c("name", "value", "metric", "dist", "params", "blockCount", "blockSizes", "blockStarts")]
  data.frame(formatted.results)
}
