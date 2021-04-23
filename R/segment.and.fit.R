#' Segment and fit a gene
#'
#' @param GENE
#'
#' @param PARAMETERS
#' @param ANNOTATION
#' @param PEAKS
#'
#' @import extraDistr
#' @export
segment.and.fit = function(
  GENE,
  PARAMETERS,
  ANNOTATION,
  PEAKS
){

  # If the gene doesn't have peaks
  if(!GENE %in% PEAKS$name){
    warn.message = paste0("No Peaks are Found for ", GENE, " in PEAKS!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(PARAMETERS))
  }

  # Use the parameters if they are defined and default to FALSE if not defined
  # plot.merged.peak <- PARAMETERS$PLOT.MERGED.PEAKS  %||% FALSE
  plot.diagnostic <- PARAMETERS$DIAGNOSTIC %||% FALSE
  write.output <- PARAMETERS$WRITE.OUTPUT %||% FALSE
  # Optional normal mixture model via mixtools
  fit.norm_mixture <- PARAMETERS$FIT.MIXTURE %||% FALSE

  # PEAKSGR
  PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = F)

  # Get Gene Information
  GENEINFO = .get.gene.anno(GENE, ANNOTATION)

  # Converting to RNA
  GENEPEAKSGR = GenomicRanges::shift(PEAKSGR, -1*GENEINFO$left+1)
  GenomicRanges::start(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::start(GENEPEAKSGR)+1]
  GenomicRanges::end(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::end(GENEPEAKSGR)]

  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = S4Vectors::split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(GenomicRanges::reduce(GENEPEAKSGR))

  # Peak Coverage
  GENEGR = GenomicRanges::GRanges(seqnames = GENEINFO$chr, IRanges::IRanges(start = 1, end = GENEINFO$exome_length), strand = GENEINFO$strand)
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(GENEGR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)

  # Segmenting & Determining which segments are peaks
  smooth.coverage = smooth.spline( BIN.COUNTS$start, BIN.COUNTS$Coverage, spar = 0.3)
  p = find.peaks(x = -smooth.coverage$y, m = 150, diff.threshold = 10^-7)

  # Formatting
  SEG.GR = generate.peaks.from.split.points(p, GENEPEAKSGR, GENEINFO, m = 100)

  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(GENEPEAKSGR, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)

  # Fitting different models
  results = data.frame()
  models = list()
  for(i in 1:length(SEG.GR)){
    seg.start = GenomicRanges::start(SEG.GR)[i]
    seg.end = GenomicRanges::end(SEG.GR)[i]
    x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]
    # Start the counts at 0 to fit the distributions better
    x.adj <- (x - seg.start) + 1e-10
    sd.scale <- sd(x.adj)
    x.scale <- x.adj / sd.scale

    # While loop
    # My thoughts are here, we can do a while loop (e.g, while end residuals > abs. residual threshold & we haven't surpassed the edge threshold),
    # continue to shrink the edges. Here's where we also place a call to the fit.residuals function
    x.range = seg.start:seg.end
    x.range.adj <- (x - seg.start) + 1e-10
    # The bin size is now 1/sd_scale
    x.range.scale <- x.range.adj / sd.scale
    # Multiple scale factor by this new bin size
    scalefactor <- length(x.scale) / sd.scale

    mod = fit.continuous.distributions(x = x.scale, sd.scale = sd.scale, seg.start = seg.start, seg.end = seg.end, fit.normal.mixture = T, max.iterations = 500)
    fits = extract.distribution.parameters(mod = mod, x = x.scale, scalefactor = scalefactor)
    fits$i = i
    # The loop ends here.

    # Adding the results to the table
    res.final = fits[fits$aic == min(fits$aic),]
    results = rbind(results, res.final)
    mod.final = mod[[res.final$dist]]
    models[[i]] = mod.final
  }

  # Acquiring a Density Distribution
  comput.fti <- function(i) { # models, SEG.GR
    # Initializing
    seg.start <- GenomicRanges::start(SEG.GR)[i]
    seg.end <- GenomicRanges::end(SEG.GR)[i]
    x = seg.start:seg.end
    x.adj <- (x - seg.start) + 1e-10
    # Extracting Model
    fti = models[[i]]
    # The bin size is now 1/sd_scale
    x.scale <- x.adj / fti$sd_scale

    if(class(fti) == "mixEM") {
      fit.data <- fti$x
      scalefactor <- length(fit.data) / fti$sd_scale
      distname <- "norm_mixture"
      dens <- dnorm_mixture(x.scale, fti)
    } else {
      fit.data <- fti$data
      # Adjust the scale of the segment
      # Multiple scale factor by this new bin size
      scalefactor <- length(fit.data) / fti$sd_scale

      # Generating data
      params <- c(as.list(fti$estimate), as.list(fti$fix.arg))
      distname <- fti$distname
      ddistname <- paste0("d", distname)
      call.params <- c(list(x = x.scale), as.list(params))
      dens <- do.call(ddistname, call.params)
    }
    dens.scale <- dens * scalefactor
    data.frame(
      "x" = x,
      "dens" = dens.scale,
      "col" = distname,
      row.names = NULL)
  }

  # Making a Nice Figure
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS) {

    distr.plotting.data = lapply(1:length(SEG.GR), comput.fti)
    # ggplot.plot(PARAMETERS=PARAMETERS, distr.plotting.data = distr.plotting.data, geneinfo=GENEINFO, bin.counts=BIN.COUNTS, seg.gr=SEG.GR, p=p)
    bpg.plot(PARAMETERS=PARAMETERS, distr.plotting.data = distr.plotting.data, geneinfo=GENEINFO, bin.counts=BIN.COUNTS, seg.gr=SEG.GR, p=p)

  }

  # Generating Peaks
  merged.peaks = SEG.GR
  S4Vectors::mcols(merged.peaks)$i = results$i
  S4Vectors::mcols(merged.peaks)$dist = results$dist
  S4Vectors::mcols(merged.peaks)$params = results$params
  S4Vectors::mcols(merged.peaks)$name = GENEINFO$gene
  merged.peaks = .rna.peaks.to.genome(merged.peaks, GENEINFO)
  GenomicRanges::start(merged.peaks) = GenomicRanges::start(merged.peaks)-1

  # Generating BED12 File
  PEAKS.FINAL = .bed6tobed12(MERGED.PEAKS = merged.peaks, ID.COLS = c("name", "i", "dist"))
  # P-Value Table
  SAMPLE.PVAL = .merge.p(PEAKSGR, MERGED.PEAKS = merged.peaks, ANNOTATION, PARAMETERS, ID.COLS = c("name", "i", "dist"))
  # Output Table
  OUTPUT.TABLE = merge(PEAKS.FINAL, SAMPLE.PVAL, by = "peak", all = T)

  return(OUTPUT.TABLE)

}
