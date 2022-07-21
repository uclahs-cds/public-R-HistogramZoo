#' Produces a GRangesList out of a GTF file, each element represents the exons of a gene or transcript
#'
#' @param gtf Path to GTF file
#' @param gene_or_transcript Whether histograms should be computed on gene annotations or transcript annotations.
#' @param select_strand Select elements belonging to a specific strand.
#' @param select_chrs Select elements on specific chromosomes.
#' @param select_ids Select elements by matching ids to genes or transcripts (depending on gene_or_transcript)
#'
#' @return A GRangesList object where each element is a GRanges object containing the exons of a gene or transcript
#'
#' @examples \dontrun{
#' gtf = system.file("extdata", "genes.gtf", package = "ConsensusPeaks")
#' gtf_gr = GTF_to_GRangesList(gtf)
#' }
#' @export
GTF_to_GRangesList = function(
  gtf,
  gene_or_transcript = c("gene", "transcript"),
  select_strand = c("*", "+", "-"),
  select_chrs = NULL,
  select_ids = NULL
){

  gene_or_transcript = match.arg(gene_or_transcript)
  select_strand = match.arg(select_strand)

  # Import GTF file
  gtf = valr::read_gtf(gtf, zero_based = F)

  # Filtering chromosomes, genes, and strand
  gtf = gtf[gtf$type == "exon",]
  if(!is.null(select_chrs)){
    gtf = gtf[gtf$chrom %in% select_chrs,]
  }
  if(select_strand != "*"){
    gtf = gtf[gtf$strand %in% select_strand,]
  }

  # Nominating the subsetting ID
  if(gene_or_transcript == "gene"){
    gtf[,"id"] = gtf[,"gene_id"]
  } else {
    gtf[,"id"] = gtf[,"transcript_id"]
  }

  # Filtering by ID
  if(!is.null(select_ids)){
    gtf = gtf[gtf$id %in% select_ids,]
  }

  # Creating a GRanges object
  gtf = gtf[,c("chrom", "start", "end", "strand", "id")]
  gtf_gr = GenomicRanges::makeGRangesFromDataFrame(df = gtf, keep.extra.columns = T)
  gtf_gr = S4Vectors::split(gtf_gr, gtf_gr$id)
  gtf_gr = GenomicRanges::reduce(gtf_gr)

  return(gtf_gr)
}
