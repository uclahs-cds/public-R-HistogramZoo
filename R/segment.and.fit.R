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
#' @param fit.mixtures Only if the method parameter is set to 'sf'. A character vector indicating distributions to fit.
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
  fit.mixtures,
  trim.peak.threshold,
  peak.length.threshold
){

  # If the gene doesn't have peaks
  if(!gene %in% peaks$name){
    warn.message = paste0("No Peaks are Found for ", gene, " in peaks!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(all.samples))
  }

  # Use the parameters if they are defined and default to FALSE if not defined
  plot.diagnostic <- diagnostic %||% FALSE
  fit.norm_mixture <- fit.mixtures %||% FALSE

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

  # Extracting Uniform Segments on Segments
  seg.gr.unif.correction = identify.uniform.segments(
    seg.gr = seg.gr,
    bin.counts = bin.counts,
    trim.peak.threshold = trim.peak.threshold,
    short.peak.threshold = peak.length.threshold
  )

  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(genepeaksgr, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)

  # Fitting different models
  results = data.frame()
  models = list()
  for(i in 1:length(seg.gr.unif.correction)){

    # Extracting data
    seg.start = GenomicRanges::start(seg.gr.unif.correction)[i]
    seg.end = GenomicRanges::end(seg.gr.unif.correction)[i]
    x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]
    x.adjusted <- (x - seg.start) + 1e-10
    x.range = seg.start:seg.end
    x.range.adjusted <- (x - seg.start) + 1e-10

    # Refit models
    mod = fit.continuous.distributions(
      x = x.adjusted,
      seg.start = seg.start,
      seg.end = seg.end,
      fit.mixtures = fit.mixtures,
      max.iterations = 500)

    # Extracting Fit Data
    fits = extract.distribution.parameters(
      mod = mod,
      x = x.adjusted)
    fits$i = seg.gr.unif.correction$i[i]

    # Adding the results to the table
    res.final = fits[fits$aic == min(fits$aic),]
    results = rbind(results, res.final)
    mod.final = mod[[res.final$dist]]
    models[[i]] = mod.final
  }

  # Making a Nice Figure
  if(gene %in% plot.merged.peaks) {
    distr.plotting.data = lapply(1:length(seg.gr.unif.correction), function(i){
      calculate.density(
        m = models[[i]],
        x = NULL,
        seg.start = GenomicRanges::start(seg.gr.unif.correction)[i],
        seg.end = GenomicRanges::end(seg.gr.unif.correction)[i],
        stepsize = 1,
        scale.density = T,
        return.df = T)
    })
    # ggplot.plot(
    #  output.dir = output.dir,
    #  distr.plotting.data = distr.plotting.data,
    #  geneinfo = geneinfo,
    #  bin.counts = bin.counts,
    #  seg.gr = seg.gr,
    #  p = p)
    bpg.plot(
      output.dir = output.dir,
      distr.plotting.data = distr.plotting.data,
      geneinfo=geneinfo,
      bin.counts=bin.counts,
      seg.gr=seg.gr,
      p=p,
      seg.gr.unif.correction = seg.gr.unif.correction,
      results = results)
  }

  # Formatting results table to include extra params
  results$width = GenomicRanges::width(seg.gr.unif.correction)
  results.dist = sapply(1:length(seg.gr), function(i) paste0(results$dist[results$i == i], collapse = ","))
  results.params = sapply(1:length(seg.gr), function(i) paste0(results$params[results$i == i], collapse = ","))
  results.mse = sapply(1:length(seg.gr), function(i){stats::weighted.mean(results$mse[results$i == i], results$width[results$i == i])})

  # Generating Peaks
  merged.peaks = seg.gr
  S4Vectors::mcols(merged.peaks)$i = 1:length(seg.gr)
  S4Vectors::mcols(merged.peaks)$dist = results.dist
  S4Vectors::mcols(merged.peaks)$params = results.params
  S4Vectors::mcols(merged.peaks)$mse = results.mse
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
