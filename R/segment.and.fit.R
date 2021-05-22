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
#' @param trim.peak.threshold See ConsensusPeaks
#' @param trim.peak.stepsize See ConsensusPeaks
#' @param residual.tolerance See ConsensusPeaks
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
  trim.peak.stepsize,
  residual.tolerance,
  eps
){

  # If the gene doesn't have peaks
  if(!gene %in% peaks$name){
    warn.message = paste0("No Peaks are Found for ", gene, " in peaks!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(all.samples))
  }

  # Use the parameters if they are defined and default to FALSE if not defined
  # plot.diagnostic <- diagnostic %||% FALSE
  # fit.norm_mixture <- fit.mixtures %||% FALSE

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
  # smooth.coverage = smooth.spline( bin.counts$start, bin.counts$Coverage, spar = 0.3)
  # p = find.peaks(x = -smooth.coverage$y, m = 150, diff.threshold = 10^-7)

  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(genepeaksgr, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)

  p = c()
  reduced.gr = GenomicRanges::reduce(genepeaksgr)
  for(i in 1:length(reduced.gr)){
    p.start = GenomicRanges::start(reduced.gr)[i]
    p.end = GenomicRanges::end(reduced.gr)[i]
    p.init = c(GenomicRanges::start(genepeaksgr), GenomicRanges::end(genepeaksgr))
    p.init = c(p.start, p.init, p.end)
    p.init = sort(unique(p.init))-p.start+1
    x = peak.counts[peak.counts >= p.start & peak.counts <= p.end]
    hist = obs.to.int.hist(x)
    p.tmp = ftc.helen(hist, p.init, eps)
    p.tmp = p.tmp+p.start-1
    p = c(p, p.tmp)
  }

  # Formatting
  seg.gr = generate.peaks.from.split.points(
    p = p,
    genepeaksgr = genepeaksgr,
    geneinfo = geneinfo,
    m = 100)

  # Fitting different models
  results = data.frame()
  models = list()
  fitted.seg.gr = GenomicRanges::GRanges()
  for(i in 1:length(seg.gr)){
    # cat(i , "\n")
    # Extracting data
    seg.start = GenomicRanges::start(seg.gr)[i]
    seg.end = GenomicRanges::end(seg.gr)[i]
    x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]
    # Adjusting X
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

    # Extract Residuals
    mod.optim = which(fits$aic == min(fits$aic))
    jc.optim.multi.distr = fits$jc[mod.optim]

    if(jc.optim.multi.distr < residual.tolerance){

      # Calculating Threshold
      peak.length = seg.end - seg.start +1
      shortest.peak =  trim.peak.threshold*peak.length
      peak.shift.max = floor((peak.length - shortest.peak)/trim.peak.stepsize)*trim.peak.stepsize

      # Initializing Refit
      refit.values = c()
      seg.start.values = seq(seg.start, seg.start + peak.shift.max, trim.peak.stepsize)
      seg.end.values = seq(seg.end - peak.shift.max, seg.end, trim.peak.stepsize)
      todo = expand.grid("seg.start.values" = seg.start.values, "seg.end.values" = seg.end.values)
      todo = todo[todo$seg.end.values - todo$seg.start.values > shortest.peak,]

      for (j in 1:nrow(todo)) {
        seg.unif.start = todo[j, "seg.start.values"]
        seg.unif.end = todo[j, "seg.end.values"]
        x.unif = peak.counts[peak.counts >=  seg.unif.start & peak.counts <= seg.unif.end]
        x.unif.adjusted = x.unif - seg.unif.start + 1e-10
        mod.unif <- fit.continuous.distributions(
          x = x.unif.adjusted,
          seg.start = seg.unif.start,
          seg.end = seg.unif.end,
          fit.mixtures = "unif",
          max.iterations = 500)
        fits.tmp = extract.distribution.parameters(
          mod = mod.unif,
          x = x.unif.adjusted)
        refit.values <- c(refit.values, fits.tmp$jc[1])
      }

      # Picking optimal fit & updating model
      max.jc = max(refit.values)
      if(max.jc > jc.optim.multi.distr){
        optim.model = which(refit.values == max.jc)
        seg.unif.start.final = todo[optim.model, "seg.start.values"][1]
        seg.unif.end.final = todo[optim.model, "seg.end.values"][1]
        x.unif.final = peak.counts[peak.counts >= seg.unif.start.final & peak.counts <= seg.unif.end.final]
        x.unif.final.adjusted = x.unif.final - seg.unif.start.final + 1e-10

        mod.final = fit.continuous.distributions(
          x = x.unif.final.adjusted,
          seg.start = seg.unif.start.final,
          seg.end = seg.unif.end.final,
          fit.mixtures = "unif",
          max.iterations = 500)
        res.final = extract.distribution.parameters(
          mod = mod.final,
          x = x.unif.final.adjusted)
        seg.gr.i = GenomicRanges::GRanges(
          seqnames = geneinfo$chr,
          IRanges::IRanges(
            start = seg.unif.start.final,
            end = seg.unif.end.final),
          strand = geneinfo$strand)
        mod.final = mod.final[[1]]
      } else {
        # Adding the results to the table
        res.final = fits[fits$aic == min(fits$aic),]
        mod.final = mod[[res.final$dist]]
        seg.gr.i = seg.gr[i]
      }
    } else {
      # Adding the results to the table
      res.final = fits[fits$aic == min(fits$aic),]
      mod.final = mod[[res.final$dist]]
      seg.gr.i = seg.gr[i]
    }
    results = rbind(results, res.final)
    models[[i]] = mod.final
    fitted.seg.gr = c(fitted.seg.gr, seg.gr.i)
  }

  # Making a Nice Figure
  if(gene %in% plot.merged.peaks) {
    distr.plotting.data = lapply(1:length(fitted.seg.gr), function(i){
      calculate.density(
        m = models[[i]],
        x = NULL,
        seg.start = GenomicRanges::start(fitted.seg.gr)[i],
        seg.end = GenomicRanges::end(fitted.seg.gr)[i],
        stepsize = 1,
        scale.density = T,
        return.df = T)
    })
    bpg.plot(
      output.dir = output.dir,
      distr.plotting.data = distr.plotting.data,
      geneinfo=geneinfo,
      bin.counts=bin.counts,
      seg.gr=seg.gr,
      p=p,
      fitted.seg.gr = fitted.seg.gr,
      results = results)
  }

  # Generating Peaks
  merged.peaks = seg.gr
  S4Vectors::mcols(merged.peaks)$i = 1:length(seg.gr)
  S4Vectors::mcols(merged.peaks)$dist = results$dist
  S4Vectors::mcols(merged.peaks)$params = results$params
  S4Vectors::mcols(merged.peaks)$mse = results$mse
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
