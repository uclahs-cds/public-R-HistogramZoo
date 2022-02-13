#' Produces a GRangesList out of a GTF file, each element represents the exons of a gene or transcript
#'
#' @param gtf.file TODO
#' @param gene.or.transcript TODO
#' @param select.strand TODO
#' @param select.chrs TODO
#' @param select.ids TODO
#'
#' @return
#' @export
gtf.to.genemodel = function(
  gtf.file,
  gene.or.transcript = c("gene", "transcript"),
  select.strand = NULL,
  select.chrs = NULL,
  select.ids = NULL
){

  gene.or.transcript = match.arg(gene.or.transcript)

  # Import GTF file
  gtf = valr::read_gtf(gtf.file, zero_based = F)

  # Filtering chromosomes, genes, and strand
  gtf = gtf[gtf$type == "exon",]
  if(!is.null(select.chrs)){
    gtf = gtf[gtf$chrom %in% select.chrs,]
  }
  if(!is.null(select.strand)){
    gtf = gtf[gtf$strand %in% select.strand,]
  }

  # Nominating the subsetting ID
  if(gene.or.transcript == "gene"){
    gtf[,"id"] = gtf[,"gene_id"]
  } else {
    gtf[,"id"] = gtf[,"transcript_id"]
  }
  if(!is.null(select.ids)){
    gtf = gtf[gtf$gene_id %in% select.ids,]
  }

  # Creating a GRanges object
  gtf = gtf[,c("chrom", "start", "end", "strand", "id")]
  gtf.gr = GenomicRanges::makeGRangesFromDataFrame(df = gtf, keep.extra.columns = T)
  gtf.gr = S4Vectors::split(gtf.gr, gtf.gr$id)
  gtf.gr = GenomicRanges::reduce(gtf.gr)

  gtf.gr
}

#' Generates a list of histogram count vectors from a coverage RLE object and a GRangesList
#'
#' @param regions TODO
#' @param coverage.rle TODO
#' @param histogram.bin.size TODO
#'
#' @return
#' @export
coverage.to.histogram = function(
  regions,
  coverage.rle,
  histogram.bin.size
){

  # Initializing
  histogram.coverage = vector("list", length(regions))
  names(histogram.coverage) = names(regions)

  # Calculating coverage
  for(i in seq_along(regions)){
    x = regions[i]
    x.name = names(x)
    x = unlist(x)
    names(x) = NULL
    bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
    bins = unlist(bins)
    GenomeInfoDb::seqlevels(bins) = GenomeInfoDb::seqlevels(coverage.rle)
    cvg = GenomicRanges::binnedAverage(
      bins = bins,
      numvar = cov,
      varname = "cvg")
    histogram.coverage[[x.name]] <- cvg$cvg
  }
  histogram.coverage
}
