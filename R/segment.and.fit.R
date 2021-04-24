#' Segment and fit a gene
#'
#' @param gene gene id from gtf file or annotation (if rna) and peaks file used to identify the set of associated peaks to merge
#' @param annotation A data frame generated using the read.gtf function for a gtf file
#' @param peaks A data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
#' @param all.samples A list of all samples with peaks in the dataset
#' @param output.tag A character string added to the names of any output files.
#' @param output.dir Output directory. If the directory does not exist, ConsensusPeaks will attempt to create the directory.
#' @param plot.merged.peaks Only if the method parameter is set to 'sf'. Either a logical value (TRUE or FALSE) indicating all or none of the merged peaks should be plotted. Otherwise, a character vector of genes whose merged peaks should be plotted.
#' @param diagnostic Only if the method parameter is set to 'sf'. A logical value indicating whether diagnostic plots for fitted distributions should be plotted.
#' @param fit.mixture Only if the method parameter is set to 'sf'. A logical value indicating whether a mixture of normal distributions should be fitted.
#'
#' @import extraDistr
#' @export
segment.and.fit = function(
  gene,
  annotation,
  peaks,
  all.samples,
  output.tag,
  output.dir,
  plot.merged.peaks,
  diagnostic,
  fit.mixture
){

  # If the gene doesn't have peaks
  if(!gene %in% peaks$name){
    warn.message = paste0("No Peaks are Found for ", gene, " in peaks!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(all.samples))
  }

  # Use the parameters if they are defined and default to FALSE if not defined
  plot.diagnostic <- diagnostic %||% FALSE
  fit.norm_mixture <- fit.mixture %||% FALSE

  # peaksgr
  peaksgr = .retrieve.peaks.as.granges(peaks = peaks, gene = gene, return.df = F)

  # Get Gene Information
  geneinfo = .get.gene.anno(gene, annotation)

  # Converting to RNA
  genepeaksgr = GenomicRanges::shift(peaksgr, -1*geneinfo$left+1)
  GenomicRanges::start(genepeaksgr) = geneinfo$DNA2RNA[GenomicRanges::start(genepeaksgr)+1]
  GenomicRanges::end(genepeaksgr) = geneinfo$DNA2RNA[GenomicRanges::end(genepeaksgr)]

  # Reduce Overlapping Peaks in the Same Sample
  genepeaksgr = S4Vectors::split(genepeaksgr, genepeaksgr$sample)
  genepeaksgr = unlist(GenomicRanges::reduce(genepeaksgr))

  # Peak Coverage
  genegr = GenomicRanges::GRanges(seqnames = geneinfo$chr, IRanges::IRanges(start = 1, end = geneinfo$exome_length), strand = geneinfo$strand)
  peak.coverage = GenomicRanges::coverage(genepeaksgr)
  bins = unlist(GenomicRanges::tile(genegr, width = 1))
  bin.counts = data.frame(GenomicRanges::binnedAverage(bins, peak.coverage, "Coverage"), stringsAsFactors = F)

  # Segmenting & Determining which segments are peaks
  smooth.coverage = smooth.spline( bin.counts$start, bin.counts$Coverage, spar = 0.3)
  p = find.peaks(x = -smooth.coverage$y, m = 150, diff.threshold = 10^-7)

  # Formatting
  seg.gr = generate.peaks.from.split.points(
    p = p,
    genepeaksgr = genepeaksgr,
    geneinfo = geneinfo,
    m = 100)

  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(genepeaksgr, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)

  # Fitting different models
  results = data.frame()
  models = list()
  for(i in 1:length(seg.gr)){
    seg.start = GenomicRanges::start(seg.gr)[i]
    seg.end = GenomicRanges::end(seg.gr)[i]
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

    mod = fit.continuous.distributions(
      x = x.scale,
      sd.scale = sd.scale,
      seg.start = seg.start,
      seg.end = seg.end,
      fit.normal.mixture = fit.mixture,
      max.iterations = 500)
    fits = extract.distribution.parameters(
      mod = mod,
      x = x.scale,
      scalefactor = scalefactor)
    fits$i = i
    # The loop ends here.

    # Adding the results to the table
    res.final = fits[fits$aic == min(fits$aic),]
    results = rbind(results, res.final)
    mod.final = mod[[res.final$dist]]
    models[[i]] = mod.final
  }

  # Acquiring a Density Distribution
  comput.fti <- function(i) { # models, seg.gr
    # Initializing
    seg.start <- GenomicRanges::start(seg.gr)[i]
    seg.end <- GenomicRanges::end(seg.gr)[i]
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
  if(gene %in% plot.merged.peaks) {
    distr.plotting.data = lapply(1:length(seg.gr), comput.fti)
    # ggplot.plot(outputdir = output.dir, distr.plotting.data = distr.plotting.data, geneinfo=geneinfo, bin.counts=bin.counts, seg.gr=seg.gr, p=p)
    bpg.plot(
      outputdir = output.dir,
      distr.plotting.data = distr.plotting.data,
      geneinfo=geneinfo,
      bin.counts=bin.counts,
      seg.gr=seg.gr,
      p=p)
  }

  # Generating Peaks
  merged.peaks = seg.gr
  S4Vectors::mcols(merged.peaks)$i = results$i
  S4Vectors::mcols(merged.peaks)$dist = results$dist
  S4Vectors::mcols(merged.peaks)$params = results$params
  S4Vectors::mcols(merged.peaks)$name = geneinfo$gene
  merged.peaks = .rna.peaks.to.genome(merged.peaks, geneinfo)
  GenomicRanges::start(merged.peaks) = GenomicRanges::start(merged.peaks)-1

  # Generating BED12 File
  peaks.final = .bed6tobed12(merged.peaks = merged.peaks, id.cols = c("name", "i", "dist"))
  # P-Value Table
  sample.pval = .merge.p(peaksgr, merged.peaks = merged.peaks, annotation, all.samples, id.cols = c("name", "i", "dist"))
  # Output Table
  output.table = merge(peaks.final, sample.pval, by = "peak", all = T)

  return(output.table)

}
