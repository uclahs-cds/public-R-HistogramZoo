find_stepfunction_chgpts = function(x){
  change_points <- which(diff(x) != 0)
  change_points_plus <- change_points+1
  keep <- (x[change_points] < x[change_points_plus])
  return( c(change_points[keep], change_points_plus[!keep]) )
}

#' Returns the indices for consecutive elements of a vector that are greater than a specified threshold
#'
#' @param x A numeric vector
#' @param threshold numeric threshold
#'
#' @return A list of coordinates with `start` and `end` coordinates
#' @export
#'
#' @examples
#' find_consecutive_threshold(c(0,0,0,1,1,1,0,0,0,1,1,1,0,0))
#' find_consecutive_threshold(c(0,0,1,2,2,0,1,1,1,0,0), threshold = 1)
find_consecutive_threshold = function(
  x,
  threshold = 0){
  x_thresholded <- rle(x > threshold)
  end_coords <- cumsum(x_thresholded$lengths)
  start_coords <- end_coords - x_thresholded$lengths + 1
  start_coords_thresholded <- start_coords[x_thresholded$values]
  end_coords_thresholded <- end_coords[x_thresholded$values]

  return(list(start = start_coords_thresholded, end = end_coords_thresholded))
}

#' segment_and_fit
#'
#' @param histogram_obj a Histogram or HistogramList object
#' @param histogram_count_threshold TODO
#' @param eps numeric hyperparameter to finetune segmentation. See `Delon et al, 2005`
#' @param seed numeric seed
#' @param truncated_models logical, whether to fit truncated distributions
#' @param max_uniform logical, whether to find a subsegment maximizing the fit of a uniform distribution
#' @param uniform_peak_threshold numeric, indicating the minimum proportion of the subsegment
#' @param uniform_peak_stepsize integer, indicating the stepsize (relative to the histogram bins) to take in the search for the uniform subsegment
#' @param remove_low_entropy logical, indicating whether to filter out low entropy regions
#' @param min_gap_size integer, indicating the minimum gap size to be filtered
#' @param histogram_metric a subset of `jaccard`, `intersection`, `ks`, `mse`, `chisq` indicating metrics to use for histogram optimization
#' @param min_peak_size integer, indication the minimum segment size
#' @param distributions a subset of `norm`, `gamma` and `unif` indicating distributions to fit
#'
#' @return A HistogramFit object representing the Histogram and results of the fit
#' @export
#' @examples \dontrun{
#' x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
#' res = segment_and_fit(x)
#' }
segment_and_fit <- function(
  histogram_obj,
  histogram_count_threshold = 0,
  eps = 1,
  seed = NULL,
  truncated_models = FALSE,
  uniform_peak_threshold = 0.75,
  uniform_peak_stepsize = 5,
  remove_low_entropy = T,
  min_gap_size = 2,
  min_peak_size = 2,
  max_uniform = T,
  histogram_metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
  distributions = c("norm", "gamma", "unif")
  ) {

  # Checking types
  stopifnot(inherits(histogram_obj, "Histogram"))
  histogram_metric <- match.arg(histogram_metric, several.ok = T)
  distributions <- match.arg(distributions, several.ok = T)

  # Extracting data
  x <- histogram_obj$histogram_data

  # Change points
  chgpts <- find_stepfunction_chgpts(x)
  # Looking for regions that surpass a hard count threshold
  x.segs <- as.data.frame(find_consecutive_threshold(x, threshold = histogram_count_threshold))
  x.segs <- x.segs[x.segs$start != x.segs$end,]

  all.points <- apply(x.segs, 1, function(segs) {
    p.init <- unname(c(segs['start'], chgpts[chgpts > segs['start'] & chgpts < segs['end']], segs['end']))
    p.init <- sort(unique(p.init)) # meaningful gaps local also needs p.init to be sorted so temporarily adding this back
    p <- ftc(x, p.init, eps)

    # Max Gap
    if(remove_low_entropy) {
      mgaps <-  meaningful_gaps_local(x = x, seg.points = p, change.points = p.init, min.gap = min_gap_size)
      p <- p[(abs(p - segs['start']) > min_peak_size & abs(p - segs['end']) > min_peak_size) | p %in% segs]
      p.pairs <- remove_max_gaps(p = p, max.gaps = mgaps, remove.short.segment = min_peak_size) # remove.short.segment can also be used to filter min_peak_size, but doesn't extend to non remove low entropy cases
    } else {
      p <- p[(abs(p - segs['start']) > min_gap_size & abs(p - segs['end']) > min_peak_size) | p %in% segs]
      p.pairs <- index_to_start_end(p)
    }

    p.pairs
  })

  # Combine the results from each segment
  all.points <- do.call('rbind.data.frame', all.points)
  all.points <- all.points[(all.points$end - all.points$start + 1) > min_peak_size,, drop = FALSE] # Current solution for min_peak_size, open to alternatives
  rownames(all.points) <- NULL # use reset.rownames?

  # Fitting different models
  models <- list()
  set.seed(seed)
  for(i in seq_len(nrow(all.points))) {
    seg <- all.points[i, ]
    seg.start <- seg[['start']]
    seg.end <- seg[['end']]
    seg.len <- seg.end - seg.start + 1
    bin.data <- x[seg.start:seg.end]

    dist.optim <- fit_distributions(bin.data, metric = histogram_metric, truncated = truncated_models, distr = distributions)
    dist.optim <- lapply(dist.optim, function(y) {
      y$seg.start <- seg.start
      y$seg.end <- seg.end
      y
    })

    # Find the maximum uniform segment
    if(max_uniform & seg.len > uniform_peak_stepsize & seg.len > ceiling(uniform_peak_threshold*seg.len)){
      max.unif.results <- identify_uniform_segment(bin.data, metric = histogram_metric, threshold = uniform_peak_threshold, step.size = uniform_peak_stepsize, max.sd.size = 0)
      # Use the maximum segment
      unif.segment <- unlist(max.unif.results[c('a', 'b')])
      bin.data.subset <- bin.data[unif.segment[1]:unif.segment[2]]
      # Fit uniform distribution on maximum uniform segment
      dist.optim.subset <- fit_distributions(bin.data.subset, metric = histogram_metric, truncated = FALSE, distr = "unif")
      # Adjust the segment starts from the shifted max uniform segment
      dist.optim.subset <- lapply(dist.optim.subset, function(y) {
        y$seg.start <- unif.segment[1] + seg.start - 1
        y$seg.end <- unif.segment[2] + seg.start - 1
        y$dist <- "unif"
        y
      })
      for(munif in names(dist.optim.subset)){dist.optim[[munif]] <- dist.optim.subset[[munif]]}
    }

    # Metric Voting
    metric.idx <- sapply(dist.optim, `[[`, "metric")
    best.models <- lapply(histogram_metric, function(met){
      metric_histograms <- dist.optim[which(metric.idx == met)]
      metric_fit <- lapply(metric_histograms, `[[`, "value")
      return( metric_histograms[which.min(metric_fit)] )
    })
    best.models <- unlist(best.models, recursive = F)
    dist.extract <- sapply(best.models, `[[`, "dist")
    dist.vote <- sort(table(dist.extract), decreasing = T)
    if(length(dist.vote) > 1 & dist.vote[1] == dist.vote[2]){
       mets <- sapply(best.models, `[[`, "metric")
       dist.majority <- dist.extract[which(mets == "jaccard")]
    } else {
      dist.majority <- names(dist.vote)[1]
    }
    best.models[['consensus']] <- dist.optim[[paste0("jaccard.", dist.majority)]]

    # Correcting for optimization via finding the minimum
    best.models <- lapply(best.models, function(mod){
      if(mod$metric %in% c("jaccard", "intersection")){ mod$value = 1 - mod$value }
      mod
    })

    best.models[['consensus']]$metric <- "consensus"
    names(best.models) <- gsub("[.].*", "", names(best.models))


    models[[i]] <- best.models
  }

  # Creating a HistogramFit object
  res <- list("models" = models, "p" = all.points,  "histogram_count_threshold" = histogram_count_threshold,
             "eps" =  eps, "seed" = seed, "truncated_models" = truncated_models, "uniform_peak_threshold" = uniform_peak_threshold,
             "uniform_peak_stepsize" = uniform_peak_stepsize, "remove_low_entropy" = remove_low_entropy, "min_gap_size" = min_gap_size,
             "min_peak_size" = min_peak_size, "max_uniform" = max_uniform, "histogram_metric" = histogram_metric, "distributions" = distributions)
  res <- c(histogram_obj, res)
  class(res) <- c("HistogramFit", class(histogram_obj))

  return(res)
}
