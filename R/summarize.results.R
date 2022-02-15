
extract.stats = function(models, gene = NULL, peak.id){
  mod = models$majority.vote
  df = data.frame(
    "peak.id" = peak.id,
    "start" = mod$seg.start,
    "end" = mod$seg.end,
    "value" = mod$value,
    "metric" = mod$metric,
    "dist" = mod$dist,
    "params" = dput.str(mod$par)
  )
  if(!is.null(gene)) {df = cbind("gene" = gene, df)}
  df
}

#' Formats results of segment.fit.agnostic
#'
#' @param result TODO
#'
#' @return
#' @export
#'
#' @examples
summarize.results.agnostic = function(
  result = NULL,
){
  res = result[['models']]
  results.tbl = lapply(seq_along(res), function(i){
    extract.stats(
      models = res[[i]],
      gene = NULL,
      peak.id = i)
  })
  results.tbl = do.call('rbind.data.frame', results.tbl)
  results.tbl
}

#' Formats results of bulk.segment.fit
#'
#' @param result TODO
#' @param output.format TODO
#'
#' TODO: Extract the correct column names for BED12
#'
#' @return
#' @export
summarize.results.bulk = function(
  result = NULL,
  output.format = c("stats.only", "bed")
){

  # Output format
  output.format = match.arg(output.format)

  # Recovering output stats
  res = result[['results']]
  res = lapply(res, `[[`, "models")
  results.tbl = lapply(seq_along(res), function(gene){
    res.gene = res[[gene]]
    lapply(seq_along(res.gene), function(i){
      extract.stats(
        models = res.gene[[i]],
        gene = names(res)[[gene]],
        peak.id = i)
    })
  })
  results.tbl = do.call('rbind.data.frame', unlist(results.tbl, recursive = F))

  if(is.null(result$gene.model) & output.format == "bed"){
    stop("A gene model is required for BED output.")
  }

  # Recovering coordinates for genes
  if(output.format == "bed"){
    gene.model = result$gene.model
    histogram.bin.size = result$histogram.bin.size
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
    coords.tbl = subset(coords.tbl, select=-c(width))
    colnames(coords.tbl)[colnames(coords.tbl) == "seqnames"] <- "chr"
    results.tbl = subset(results.tbl, select=-c(start, end))
    results.tbl = cbind.data.frame(results.tbl, coords.tbl)
  }

  results.tbl
}
