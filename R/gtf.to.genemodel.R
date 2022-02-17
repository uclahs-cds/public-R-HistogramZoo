#' Produces a GRangesList out of a GTF file, each element represents the exons of a gene or transcript
#'
#' @param gtf.file Path to GTF file
#' @param gene.or.transcript Whether histograms should be computed on gene annotations or transcript annotations.
#' @param select.strand Select elements belonging to a specific strand.
#' @param select.chrs Select elements on specific chromosomes.
#' @param select.ids Select elements by matching ids to genes or transcripts (depending on gene.or.transcript)
#'
#' @return A GRangesList object where each element is a GRanges object containing the exons of a gene or transcript
#'
#' @examples \dontrun{
#' gtf.file = system.file("extdata", "genes.gtf", package = "ConsensusPeaks")
#' gtf.gr = gtf.to.genemodel(gtf.file)
#' }
#' @export
gtf.to.genemodel = function(
  gtf.file,
  gene.or.transcript = c("gene", "transcript"),
  select.strand = c(".", "+", "-"),
  select.chrs = NULL,
  select.ids = NULL
){

  gene.or.transcript = match.arg(gene.or.transcript)
  select.strand = match.arg(select.strand)

  # Import GTF file
  gtf = valr::read_gtf(gtf.file, zero_based = F)

  # Filtering chromosomes, genes, and strand
  gtf = gtf[gtf$type == "exon",]
  if(!is.null(select.chrs)){
    gtf = gtf[gtf$chrom %in% select.chrs,]
  }
  if(select.strand != "."){
    gtf = gtf[gtf$strand %in% select.strand,]
  }

  # Nominating the subsetting ID
  if(gene.or.transcript == "gene"){
    gtf[,"id"] = gtf[,"gene_id"]
  } else {
    gtf[,"id"] = gtf[,"transcript_id"]
  }

  # Filtering by ID
  if(!is.null(select.ids)){
    gtf = gtf[gtf$id %in% select.ids,]
  }

  # Creating a GRanges object
  gtf = gtf[,c("chrom", "start", "end", "strand", "id")]
  gtf.gr = GenomicRanges::makeGRangesFromDataFrame(df = gtf, keep.extra.columns = T)
  gtf.gr = S4Vectors::split(gtf.gr, gtf.gr$id)
  gtf.gr = GenomicRanges::reduce(gtf.gr)

  gtf.gr
}
