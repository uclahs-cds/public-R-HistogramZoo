
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
  split.by.strand = TRUE,
  select.chrs = NULL,
  select.genes = NULL
){

  # Import GTF file
  gtf = valr::read_gtf(gtf.file, zero_based = F)

  # Filtering exons, chromosomes, genes
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

  # Creating a GRanges object
  seq.info = if(is.null(select.chrs)) sort(unique(gtf$chrom)) else select.chrs
  gtf = gtf[,c("chrom", "start", "end", "strand", "id")]
  gtf.gr = GenomicRanges::makeGRangesFromDataFrame(
    df = gtf,
    keep.extra.columns = T,
    seqinfo = seq.info
  )

  # Splitting by strand
  if(split.by.strand){
    gtf.gr = split.strand(gr = gtf.gr)
    gtf.gr = split.id(gr = gtf.gr, reduce = T)
  } else {
    gtf.gr = split.id(gr = gtf.gr, reduce = T)
  }

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
#' @param split.by.strand
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
#' TODO: Incorporate strand splitting widely
bed.to.hist = function(
  filenames,
  n_fields = 3,
  regions.of.interest = NULL,
  split.by.strand = TRUE,
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  histogram.bin.size = 1,
  ...
){

  # Load BED files & convert to base 1
  segs = valr::read_bed( filenames, n_fields = n_fields, ...)
  if(n_fields == 12){ segs = valr::bed12_to_exons( segs ) }
  segs.gr = GenomicRanges::makeGRangesFromDataFrame( segs )
  segs.gr = base0.to.base1(segs.gr)

  # Creating a GRanges list
  if(!is.null(regions.of.interest)){
    if( is(regions.of.interest, "GRanges") ){
      # Selecting specific genomic regions
      regions.gr = GenomicRanges::split(regions.of.interest, seq_along(regions.of.interest))
    } else if (is.character(regions.of.interest) & ! is.null(gtf.file) ){
      # RNA coordinates
      regions.gr = gtf.to.genemodel(
        gtf.file = gtf.file,
        split.by.strand = split.by.strand,
        gene.or.transcript = gene.or.transcript,
        select.chrs = levels(GenomicRanges::seqnames(segs.gr)),
        select.genes = regions.of.interest)
      # Filter regions (maybe)
    }
  } else if(!is.null(gtf.file)){
    # RNA coordinates
    regions.gr = gtf.to.genemodel(
      gtf.file = gtf.file,
      split.by.strand = split.by.strand,
      gene.or.transcript = gene.or.transcript,
      select.chrs = levels(GenomicRanges::seqnames(segs.gr)))
    # Filter regions (maybe)
  } else {
    # Genome coordinates
    regions.gr = GenomicRanges::reduce(segs.gr)
    regions.gr = GenomicRanges::split(regions.gr, seq_along(regions.gr))
  }

  if(split.by.strand){
    segs.strand = split.strand(segs.gr)
    segs.coverage = lapply(segs.strand, GenomicRanges::coverage)
    common.strand = intersect(names(segs.coverage), names(regions.gr))
    histogram.coverage = lapply(common.strand, function(st){
      coverage.to.hist(
        cov = segs.coverage[[st]],
        regions = regions.gr[[st]],
        histogram.bin.size = histogram.bin.size
      )
    })
    histogram.coverage = unlist(histogram.coverage, recursive = F)
    regions.gr = do.call(c, lapply(common.strand, function(i) regions.gr[[i]]))
  } else {
    segs.coverage = GenomicRanges::coverage( segs.gr )
    histogram.coverage = coverage.to.hist(
      cov = segs.coverage,
      regions = regions.gr,
      histogram.bin.size = histogram.bin.size)
  }

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
