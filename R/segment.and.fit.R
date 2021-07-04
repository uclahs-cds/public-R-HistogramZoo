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
  truncated.models = FALSE,
  uniform.peak.threshold = 0.75,
  uniform.peak.stepsize = 5,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
  eps = 1
){
  histogram.metric = match.arg(histogram.metric)

  # If the gene doesn't have peaks
  if(!gene %in% peaks$name){
    warn.message = paste0("No Peaks are Found for ", gene, " in peaks!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(all.samples))
  }

  # Use the parameters if they are defined and default to FALSE if not defined
  # plot.diagnostic = diagnostic %||% FALSE
  # fit.norm_mixture = fit.mixtures %||% FALSE

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
    # cat(i, "\n")
    p.start = GenomicRanges::start(reduced.gr)[i]
    p.end = GenomicRanges::end(reduced.gr)[i]
    tmp.gr = genepeaksgr[S4Vectors::subjectHits(GenomicRanges::findOverlaps(reduced.gr[i], genepeaksgr))]
    p.init = c(GenomicRanges::start(tmp.gr), GenomicRanges::end(tmp.gr))
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
    seg.len = seg.end - seg.start + 1
    x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]
    # Adjusting X
    x.adjusted = (x - seg.start) + 1
    #x.range = seg.start:seg.end
    #x.range.adjusted = (x - seg.start) + 1e-10

    bin.data = obs.to.int.hist(x.adjusted, as.df = TRUE, add.zero.endpoints = FALSE)
    dist.optim = fit.distributions.optim(bin.data, metric = histogram.metric, truncated = truncated.models)
    dist.optim = lapply(dist.optim, function(y) {
      y$seg.start = seg.start
      y$seg.end = seg.end
      y
    })

    # Find the maximum uniform segment
    if(seg.len > uniform.peak.stepsize & seg.len > ceiling(uniform.peak.threshold*seg.len)){
      max.unif.results = find.uniform.segment(bin.data$Freq, threshold = uniform.peak.threshold, step.size = uniform.peak.stepsize, max.sd.size = 0)
      # Use the maximum segment
      unif.segment = unlist(max.unif.results[c('a', 'b')])
      unif.segment.adj =  unif.segment
      x.subset = x[x >= unif.segment.adj[1] & x <= unif.segment.adj[2]]
      bin.data.subset = bin.data[bin.data$x >= unif.segment.adj[1] & bin.data$x <= unif.segment.adj[2],]
      bin.data.subset$x = bin.data.subset$x - unif.segment.adj[1] + 1
      # Fit uniform distribution on maximum uniform segment
      dist.optim.subset = fit.distributions.optim(bin.data.subset, metric = histogram.metric, truncated = FALSE, distr = "unif")
      # Adjust the segment starts from the shifted max uniform segment
      dist.optim.subset$unif$seg.start = unif.segment.adj[1] + seg.start
      dist.optim.subset$unif$seg.end = unif.segment.adj[2] + seg.start

      dist.optim$max_unif = dist.optim.subset$unif
    }

    # Extract all of the "value" keys from each of the models
    # The value = the metric we were optimizing
    metric.values = unlist(lapply(dist.optim, `[[`, "value"))

    mod.final = dist.optim[[which.min(metric.values)]]

    seg.gr.i = GenomicRanges::GRanges(
             seqnames = geneinfo$chr,
             IRanges::IRanges(
               start = mod.final$seg.start,
               end = mod.final$seg.end),
             strand = geneinfo$strand)

    models[[i]] = mod.final
    fitted.seg.gr = c(fitted.seg.gr, seg.gr.i)
  }

  results = do.call(rbind.data.frame, lapply(models, function(m) {
     if(histogram.metric %in% c("jaccard", "intersect")) {
       metric = 1 - m$value
     } else {
       metric = m$value
     }
    data.frame(dist = m$dist, params = dput.str(m$par), metric = metric)
  }))


  # Making a Nice Figure
  if(gene %in% plot.merged.peaks) {
    distr.plotting.data = lapply(models, function(m) {
      x = seq(m$seg.start + 1, m$seg.end, by = 1)
      # Since the density is from 1 to max of segment, pass in the sequence 1...max
      dens = m$dens(seq_along(x) + 1)
      data.frame(
        x = x,
        dens = dens,
        dist = m$dist
      )
    })

    bpg.plot(
      output.tag = output.tag,
      output.dir = output.dir,
      distr.plotting.data = distr.plotting.data,
      geneinfo=geneinfo,
      bin.counts=bin.counts,
      seg.gr=seg.gr,
      p=p,
      fitted.seg.gr = fitted.seg.gr,
      results = results,
      histogram.metric = histogram.metric)
  }

  # Generating Peaks
  merged.peaks = seg.gr
  S4Vectors::mcols(merged.peaks)$i = 1:length(seg.gr)
  S4Vectors::mcols(merged.peaks)$dist = results$dist
  S4Vectors::mcols(merged.peaks)$params = results$params
  S4Vectors::mcols(merged.peaks)$metric = results$metric
  S4Vectors::mcols(merged.peaks)$name = geneinfo$gene
  merged.peaks = .rna.peaks.to.genome(merged.peaks, geneinfo)
  GenomicRanges::start(merged.peaks) = GenomicRanges::start(merged.peaks)-1

  # Generating BED12 File
  peaks.final = .bed6tobed12(
    merged.peaks = merged.peaks,
    name.id = "name",
    peak.id = "i",
    score.id = "metric",
    distribution.id = "dist",
    meta.id = "params"
  )

  # P-Value Table
  sample.pval = .merge.p(
    peaksgr = peaksgr,
    merged.peaks = merged.peaks,
    annotation = annotation,
    all.samples = all.samples,
    name.id = "name",
    peak.id = "i"
  )

  # Output Table
  output.table = merge(peaks.final, sample.pval, by = "peak", all = T)
  return(output.table)
  # return(list(output.table = output.table, results = results))

}
