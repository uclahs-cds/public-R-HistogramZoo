#' Generates GenomicHistogram objects from the coverage of a bigWig file 
#'
#' @param filename Name of bigwig file for import
#' @param strand The strand of bigwig file from which the data originates. Default "*". If strand is "+" or "-", the strand will also be used to select regions of matching strand.
#' @param score_threshold A hard threshold for the score of the bigwig file. Scores higher than the threshold will be used in the computation of the histogram.
#' @param regions A GRanges or GRangesList object representing regions of interest defining histograms, a GRangesList object allows for the specification of non-continuous segments,
#' default NULL - regions are taken as the set of elements with nonzero coverage
#' @param gtf A GTF file to select regions of interest. Default NULL
#' @param histogram_bin_size The bin size (base-pairs) to bin signal into a histogram. Default 50
#' @param gene_or_transcript Whether histograms should be computed on gene annotations or transcript annotations. Default gene
#' @param select_chrs Select elements on specific chromosomes. Default NULL
#' @param select_ids Select elements by matching ids to genes or transcripts (depending on gene_or_transcript). Default NULL
#'
#' @return A list of GenomicHistogram objects
#'
#' @examples \dontrun{
#' filename = system.file("extdata", "bigwigs",  "S1.bw", package = "ConsensusPeaks")
#' regions = GenomicRanges::GRanges(
#' seqnames = "chr1",
#' IRanges::IRanges(start = c(17950, 19350),
#'                  end = c(18000, 19600)),
#' strand = "*")
#'
#' histograms = bigWig_to_histogram(
#' filename = filename,
#' regions = regions,
#' histogram_bin_size = 10)
#' }
#'
#'
#' @export
bigWig_to_histogram = function(
  filename,
  strand  = c("*", "+", "-"),
  score_threshold = 0,
  regions = NULL,
  gtf = NULL,
  histogram_bin_size = 50,
  gene_or_transcript = c("gene", "transcript"),
  select_chrs = NULL,
  select_ids = NULL
){
  
  # Error checking
  if(!file.exists(filename) | length(filename) > 1){
    stop(filename, " doesn't exist")
  }
  if(length(filename) > 1){
    stop("filename needs to be length 1")
  }
  if(!is.numeric(score_threshold)){
    stop("score threshold must be numeric")
  }
  if(!(inherits(regions, "GRanges") | inherits(regions, "CompressedGRangesList") | is.null(regions))){
    stop("regions must be a GRanges or GRangesList object or set to NULL")
  }
  if(!is_equal_integer(histogram_bin_size) | histogram_bin_size < 1 ){
    stop("histogram_bin_size must be a positive integer")
  }
  if(is.null(gtf) & is.null(regions)){
    warning("computing regions on coverage data, this may take a long a time.")
  }

  # Load BigWig file
  bigwig = valr::read_bigwig(filename, set_strand = strand)
  bigwig_gr = GenomicRanges::makeGRangesFromDataFrame(bigwig, keep.extra.columns = T)
  bigwig_gr = bigwig_gr[bigwig_gr$score > score_threshold]
  bigwig_coverage = GenomicRanges::coverage(bigwig_gr, weight = "score")

  # Loading regions
  if(inherits(regions, "GRanges")){
    # Filtering by strand
    if(strand %in% c("+", "-")){
      regions = regions[GenomicRanges::strand(regions) == strand]
    }
    # Generate names for regions
    ids = if(!is.null(regions$name)) regions$name else generate_identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  } else if (inherits(regions, "CompressedGRangesList")){
    # Filtering by strand
    region_strand = as.character(
      sapply(regions, function(x) attributes(strand(x))$values)
    )
    regions = regions[region_strand == strand]
    # Generate names for regions
    if(is.null(names(regions))){
      ids <- generate_identifiers( unlist( range(regions) ) )
      names(regions) <- ids
    }
  } else if(!is.null(gtf)){
    regions = GTF_to_GRangesList(
      gtf = gtf,
      gene_or_transcript = gene_or_transcript,
      select_strand = if(strand %in% c("+", "-")) strand else NULL,
      select_chrs = select_chrs,
      select_ids = select_ids
    )
  } else { # if no user-specified regions are provided
    regions = GenomicRanges::reduce(bigwig_gr)
    ids = generate_identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  }

  # Sort regions
  regions = sort(regions)
  
  # Generating Histogram objects
  return(
    structure(
      lapply(names(regions), function(id) coverage_to_histogram(
        region = regions[[id]],
        region_id = id,
        coverage = bigwig_coverage,
        histogram_bin_size = histogram_bin_size
        )
      ),
      names = names(regions)
    )
  )
}
