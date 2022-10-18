#' Generate GenomicHistogram objects from BED files using coordinates defined by a GTF file
#'
#' @param filenames A vector of BED filenames. The `name` column of the BED files must indicate gene or transcript name
#' @param n_fields Number of columns in the BED file that conform to BED file standards, default 4
#' @param gtf A GTF file
#' @param histogram_bin_size The bin size (base-pairs) to bin signal into a histogram, default 1
#' @param gene_or_transcript `gene` or `transcript` indicating whether histograms should be built on gene coordinates or transcript coordinates, default gene
#' @param select_strand `*`, `+` or `-`  to filter regions by strand, default '*'
#' @param select_chrs a vector of chromosomes to filter regions by chromosome, default NULL
#' @param select_ids gene or transcript ids to filter regions by gene or transcript name, must correspond to GTF ids, defualt NULL
#' @param allow_overlapping_segments_per_sample logical, if FALSE, overlapping segments in the same region in the same file will be de-duplicated in the coverage calculation, default NULL
#' NOTE: regions are determined by GTF gene/transcript IDs and the name column of BED files. If TRUE, they will be taken as separate input, default FALSE
#'
#' @return a list of GenomicHistogram objects
#'
#' @examples \dontrun{
#' datadir <- system.file("extdata", "rna_bedfiles", package = "HistogramZoo")
#' filenames <- file.path(datadir, paste0("Sample.", 1:20, ".bed"))
#' gtf <- system.file("extdata", "genes.gtf", package = "HistogramZoo")
#'
#' histograms <- transcript_BED_to_histogram(
#' filenames = filenames,
#' n_fields = 12,
#' gtf = gtf,
#' gene_or_transcript = "gene",
#' histogram_bin_size = 10)
#' }
#'
#' @export
transcript_BED_to_histogram = function(
  filenames,
  n_fields = c(4, 6, 12),
  gtf,
  histogram_bin_size = 1,
  gene_or_transcript = c("gene", "transcript"),
  select_strand = c("*", "+", "-"),
  select_chrs = NULL,
  select_ids = NULL,
  allow_overlapping_segments_per_sample = F
){

  # Error checking
  existing_files <- file.exists(filenames)
  if(!all(existing_files)){
    stop(paste0(basename(filenames[!existing_files]), collapse = ","), " don't exist")
  }
  if(!is_equal_integer(histogram_bin_size) | histogram_bin_size < 1 ){
    stop("histogram_bin_size must be a positive integer")
  }
  if(!is.logical(allow_overlapping_segments_per_sample)){
    stop("allow_overlapping_segments_per_sample must be logical")
  }

  # Reading in peak files
  peaks <- lapply(filenames, function(filename){
    segs <- valr::read_bed( filename, n_fields = n_fields, comment = "#")
    if(n_fields == 12){ segs <- valr::bed12_to_exons( segs ) }
    segs_gr <- GenomicRanges::makeGRangesFromDataFrame( segs, keep.extra.columns = T )
    segs_gr <- base0_to_base1(segs_gr)
    if(!allow_overlapping_segments_per_sample){
      segs_gr <- S4Vectors::split(segs_gr, f = segs_gr$name)
      segs_gr <- GenomicRanges::reduce(segs_gr)
      segs_gr <- unlist(segs_gr)
    }
    segs_gr
  })
  peaks <- do.call(c, peaks)
  peaks <- split(peaks, f = names(peaks))

  # Loading GTF file
  regions <- GTF_to_GRangesList(
    gtf = gtf,
    gene_or_transcript = gene_or_transcript,
    select_strand = select_strand,
    select_chrs = select_chrs,
    select_ids = select_ids
  )

  ids <- intersect(names(peaks), names(regions))
  if(length(ids) == 0) {
    warning("No intersecting IDs between regions and peaks!")
    return(NULL)
  }

  # Generating Histogram objects
  return(
    structure(
      lapply(ids, function(id) coverage_to_histogram(
        region = regions[[id]],
        region_id = id,
        coverage = GenomicRanges::coverage(peaks[[id]]),
        histogram_bin_size = histogram_bin_size
        )
      ),
      names = ids
    )
  )
}
