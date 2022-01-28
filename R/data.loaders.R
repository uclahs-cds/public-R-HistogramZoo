

#' Produces a GRangesList out of a GTF file, each element represents the exons of a gene or transcript
#'
#' @param gtf.file
#' @param gene.or.transcript
#' @param select.chrs
#'
#' @return
#' @export
#'
#' @examples
gtf.to.genemodel = function(
  gtf.file,
  gene.or.transcript = c("gene", "transcript"),
  select.chrs = NULL,
  select.genes = NULL
){

  # Import GTF file
  gtf = valr::read_gtf(gtf.file)

  # Filtering
  gtf = gtf[gtf$type == "exon",]
  if(!is.null(select.chrs)){
    gtf = gtf[gtf$chrom %in% select.chrs,]
  }
  if(!is.null(select.genes)){
    gtf = gtf[gtf$gene_id %in% select.genes,]
  }

  # Nominating the subsetting ID
  if(gene.or.transcript == "gene"){
    gtf[,"id"] = gtf[,"gene_id"]
  } else {
    gtf[,"id"] = gtf[,"transcript_id"]
  }

  seq.info = if(is.null(select.chrs)) sort(unique(gtf$chrom)) else select.chrs
  gtf = gtf[,c("chrom", "start", "end", "strand", "id")]
  gtf.gr = GenomicRanges::makeGRangesFromDataFrame(
    gtf,
    keep.extra.columns = T,
    seqinfo = seq.info
  )
  gtf.gr = S4Vectors::split(gtf.gr, gtf.gr$id)
  gtf.gr = GenomicRanges::reduce(gtf.gr)
  gtf.gr
}

#' Generates a list of histogram count vectors from a coverage RLE object and a GRangesList
#'
#' @param regions
#' @param cov
#' @param histogram.bin.size
#'
#' @return
#' @export
#'
#' @examples
#' TODO: Implement error checking
coverage.to.hist = function(
  regions,
  cov,
  histogram.bin.size
){
  histogram.coverage = list()
  for(i in seq_along(regions)){
    x = regions[i]
    x.name = names(x)
    x = unlist(x)
    bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
    bins = unlist(bins)
    cvg = GenomicRanges::binnedAverage(
      bins = bins,
      numvar = cov,
      varname = "cvg")
    cvg = cvg$cvg
    names(cvg) = NULL
    histogram.coverage[[x.name]] <- cvg
  }
  histogram.coverage
}

#' Imports BED files and exports a list of count vectors, gene models and a histogram bin size
#'
#' @param filenames
#' @param n_fields
#' @param regions.of.interest
#' @param gtf.file
#' @param gene.or.transcript
#' @param histogram.bin.size
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' TODO: Add error checks
#' TODO: Allow uses to offer a vector of gene names or a GRanges object of interested regions
bed.to.hist = function(
  filenames,
  n_fields = 3,
  regions.of.interest = NULL,
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  histogram.bin.size = 1,
  ...
){

  # Load BED files
  segs = valr::read_bed( filenames, n_fields = n_fields, ...)
  if(n_fields == 12){ segs = valr::bed12_to_exons( segs ) }
  segs.gr = GenomicRanges::makeGRangesFromDataFrame( segs )

  # Calculating coverage using GenomicRanges
  segs.coverage = GenomicRanges::coverage( segs.gr )

  # Creating a GRanges list
  if(!is.null(regions.of.interest)){
    if( is(regions.of.interest, "GRanges") ){
      # Selecting specific genomic regions
      regions.gr = regions.of.interest
    } else if (is.character(regions.of.interest) & ! is.null(gtf.file) ){
      # RNA coordinates
      regions.gr = gtf.to.genemodel(
        gtf.file = gtf.file,
        gene.or.transcript = gene.or.transcript,
        select.chrs = names(segs.coverage),
        select.genes = regions.of.interest)
      select.regions = S4Vectors::subjectHits(GenomicRanges::findOverlaps(segs.gr, regions.gr))
      regions.gr = regions.gr[sort(unique(select.regions))]
    }
  } else if(!is.null(gtf.file)){
    # RNA coordinates
    regions.gr = gtf.to.genemodel(
      gtf.file = gtf.file,
      gene.or.transcript = gene.or.transcript,
      select.chrs = names(segs.coverage))
    select.regions = S4Vectors::subjectHits(GenomicRanges::findOverlaps(segs.gr, regions.gr))
    regions.gr = regions.gr[sort(unique(select.regions))]
  } else {
    # Genome coordinates
    regions.gr = GenomicRanges::reduce(segs.gr)
    regions.gr = GenomicRanges::split(regions.gr, seq_along(regions.gr))
  }

  # Generate coverage vectors
  histogram.coverage = coverage.to.hist(
    cov = segs.coverage,
    regions = regions.gr,
    histogram.bin.size = histogram.bin.size)

  # Returning list of histograms, region coordinates
  list(
    histogram.coverage = histogram.coverage,
    gene.model = regions.gr,
    histogram.bin.size = histogram.bin.size)
}

# Imports BigWig files and exports something usable by the segmentation function
bigwig.to.hist = function(
  filename,
  set.strand = c("+", "-"),
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  histogram.bin.size = 1
){

  # Load BigWig file
  bigwig = valr::read_bigwig(filename, set_strand = set.strand)
  bigwig.gr = GenomicRanges::makeGRangesFromDataFrame(bigwig, keep.extra.columns = T)
  segs.coverage = GenomicRanges::coverage(bigwig.gr, weight = "score")
  bigwig.gr = bigwig.gr[bigwig.gr$score != 0]

  # Loading regions
  if(!is.null(gtf.file)){
    regions.gr = gtf.to.genemodel(
      gtf.file = gtf.file,
      gene.or.transcript = gene.or.transcript,
      select.chrs = names(segs.coverage))

    # Filtering regions.gr
    select.regions = S4Vectors::subjectHits(GenomicRanges::findOverlaps(bigwig.gr, regions.gr))
    regions.gr = regions.gr[sort(unique(select.regions))]
  } else {
    # Genome coordinates
    # Rethink this a little bit***
    # regions.gr = GenomicRanges::reduce(cov.gr)
    # regions.gr = GenomicRanges::split(regions.gr, 1:length(regions.gr))
  }

  histogram.coverage = coverage.to.hist(
    cov = segs.coverage,
    regions = regions.gr,
    histogram.bin.size = histogram.bin.size)

  # Returning list of histograms, region coordinates
  list(
    histogram.coverage = histogram.coverage,
    gene.model = regions.gr,
    histogram.bin.size = histogram.bin.size)
}
