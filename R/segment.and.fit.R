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
  eps = 1,
  remove.low.entropy = T,
  max.uniform = T
){
  histogram.metric = match.arg(histogram.metric, several.ok = T)

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

  # Looks for the lower changepoint of a step function
  find.stepfunction.chgpts = function(x){
    chg.pts = which(diff(x) != 0)
    chg.pts.plus = chg.pts+1
    keep.chg.pts = (x[chg.pts] < x[chg.pts.plus])
    c(chg.pts[keep.chg.pts], chg.pts.plus[!keep.chg.pts])
  }

  p = c()
  max.gaps = data.frame()
  reduced.gr = GenomicRanges::reduce(genepeaksgr)
  for(i in 1:length(reduced.gr)){
    # cat(i, "\n")

    # Extracting Coverage
    p.start = GenomicRanges::start(reduced.gr)[i]
    p.end = GenomicRanges::end(reduced.gr)[i]
    bin.data = bin.counts$Coverage[p.start:p.end]

    # Change points
    p.init = find.stepfunction.chgpts(bin.data)
    p.init = sort(unique(c(1, p.init, p.end - p.start + 1)))

    # FTC
    p.tmp = ftc.helen(bin.data, p.init, eps)

    # Max Gap
    if(remove.low.entropy){
      mgaps = meaningful.gaps.local(x = bin.data, seg.points = p.tmp, change.points = p.init)
      mgaps$Var1 = mgaps$Var1+p.start-1
      mgaps$Var2 = mgaps$Var2+p.start-1
      max.gaps = rbind(max.gaps, mgaps[,c("Var1", "Var2")])
    }

    # Updating
    p.tmp = p.tmp+p.start-1
    p = c(p, p.tmp)

  }

  # Formatting
  seg.gr = generate.peaks.from.split.points(
    p = p,
    genepeaksgr = genepeaksgr,
    geneinfo = geneinfo,
    m = 100)

  # Removing Max Gaps
  if(remove.low.entropy){
    seg.gr = remove.max.gaps(
      geneinfo = geneinfo,
      seg.gr = seg.gr,
      max.gaps = max.gaps,
      remove.short.segment = 100)
  }

  # Fitting different models
  results = data.frame()
  models = list()
  for(i in 1:length(seg.gr)){
    # cat(i , "\n")
    # Extracting data
    seg.start = GenomicRanges::start(seg.gr)[i]
    seg.end = GenomicRanges::end(seg.gr)[i]
    seg.len = seg.end - seg.start + 1
    bin.data = bin.counts$Coverage[seg.start:seg.end]

    dist.optim = fit.distributions.optim(bin.data, metric = histogram.metric, truncated = truncated.models)
    dist.optim = lapply(dist.optim, function(y) {
      y$seg.start = seg.start
      y$seg.end = seg.end
      y
    })

    # Find the maximum uniform segment
    if(max.uniform & seg.len > uniform.peak.stepsize & seg.len > ceiling(uniform.peak.threshold*seg.len)){
      max.unif.results = find.uniform.segment(bin.data, metric = histogram.metric, threshold = uniform.peak.threshold, step.size = uniform.peak.stepsize, max.sd.size = 0)
      # Use the maximum segment
      unif.segment = unlist(max.unif.results[c('a', 'b')])
      bin.data.subset = bin.data[unif.segment[1]:unif.segment[2]]
      # Fit uniform distribution on maximum uniform segment
      dist.optim.subset = fit.distributions.optim(bin.data.subset, metric = histogram.metric, truncated = FALSE, distr = "unif")
      # Adjust the segment starts from the shifted max uniform segment
      dist.optim.subset = lapply(dist.optim.subset, function(y) {
        y$seg.start = unif.segment[1] + seg.start
        y$seg.end = unif.segment[2] + seg.start
        y$dist = "unif"
        y
      })
      for(munif in names(dist.optim.subset)){dist.optim[[munif]] <- dist.optim.subset[[munif]]}
    }

    # Metric Voting
    value.df = data.frame(
      "i" = i,
      "value" = unlist(lapply(dist.optim, `[[`, "value")),
      "dist" = unlist(lapply(dist.optim, `[[`, "dist")),
      "metric" = unlist(lapply(dist.optim, `[[`, "metric")),
      "params" = unlist(lapply(dist.optim, function(m) dput.str(m$par))),
      "seg.start" = unlist(lapply(dist.optim, `[[`, "seg.start")),
      "seg.end" = unlist(lapply(dist.optim, `[[`, "seg.end")),
      "final" = 0,
      stringsAsFactors=F
    )
    distr.vote = aggregate(value ~ metric, value.df, FUN = min)
    vote.df = merge(value.df, distr.vote, by = c("metric", "value"))
    distr.tally = table(vote.df$dist)
    best.distr = ifelse(sum(distr.tally == max(distr.tally))>1, vote.df$dist[vote.df$metric == "jaccard"], names(distr.tally)[which.max(distr.tally)])
    final.res = value.df[value.df$metric == "jaccard" & value.df$dist == best.distr,]
    final.res$final <- 1
    vote.df = rbind.data.frame(vote.df, final.res)

    mod.final = dist.optim[[paste0("jaccard.", best.distr)]]
    models[[i]] = mod.final
    results = rbind(results, vote.df)
  }

  # Correcting for optimization via finding the minimum
  results$value[results$metric %in% c("jaccard", "intersection")] <- 1 - results$value[results$metric %in% c("jaccard", "intersection")]

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
      p=p,
      results = results)
  }

  # Generating Peaks
  results$chr = geneinfo$chr
  results$name = geneinfo$gene
  results$strand = geneinfo$strand
  results = results[results$final == 1,,drop = F]
  merged.peaks = GenomicRanges::makeGRangesFromDataFrame(results, keep.extra.columns = T)
  merged.peaks = .rna.peaks.to.genome(merged.peaks, geneinfo)
  GenomicRanges::start(merged.peaks) = GenomicRanges::start(merged.peaks)-1

  # Generating BED12 File
  peaks.final = .bed6tobed12(
    merged.peaks = merged.peaks,
    name.id = "name",
    peak.id = "i",
    score.id = "value",
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

}
