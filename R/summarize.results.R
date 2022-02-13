
extract.stats = function(models, gene, peak.id){
  mod = models$majority.vote
  data.frame(
    "gene" = gene,
    "peak.id" = peak.id,
    "start" = mod$seg.start,
    "end" = mod$seg.end,
    "value" = mod$value,
    "metric" = mod$metric,
    "dist" = mod$dist,
    "params" = dput.str(mod$par)
  )
}

# TODO: Create a reasonable format for basic histograms
#' Formats results
#'
#' @param results.object TODO
#' @param output.format TODO
#'
#' @return
#' @export
summarize.results = function(
  results.object,
  output.format = c("stats.only", "bed")
){

  output.format = match.arg(output.format)
  if(is.null(results.object$gene.model) & output.format == "bed"){
    stop("A gene model is required for BED output.")
  }

  # Recovering output stats
  results = results.object[['results']]
  results = lapply(results, `[[`, "models")
  results.tbl = lapply(seq_along(results), function(gene){
    res.gene = results[[gene]]
    lapply(seq_along(res.gene), function(i){
      extract.stats(
        models = res.gene[[i]],
        gene = names(results)[[gene]],
        peak.id = i)
    })
  })
  results.tbl = do.call('rbind.data.frame', unlist(results.tbl, recursive = F))

  # Recovering coordinates
  if(output.format == "bed"){
    gene.model = results.object$gene.model
    histogram.bin.size = results.object$histogram.bin.size
    bed.coords = apply(results.tbl, 1, function(peak){
      gene = peak['gene']
      x = unlist(gene.model[gene])
      bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
      bins = unlist(bins)
      coords.gr = GenomicRanges::reduce(bins[peak['start']:peak['end']])
      coords.gr = base1.to.base0(coords.gr)
      coords.gr = bed6tobed12(coords.gr)
      coords.gr
    })
    coords.tbl = do.call(c, bed.coords)
    coords.tbl = data.frame(coords.tbl)
    # TODO: Extract the correct column names for BED12
    # TODO: Remove the start and end for stats
    results.tbl = cbind.data.frame(results.tbl, coords.tbl)
  }

  results.tbl
}
