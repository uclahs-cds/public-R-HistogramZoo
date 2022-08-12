#' Find changepoints in a vector with uniform stretches of values
#' @param x A numeric vector
find_stepfunction_chgpts <- function(x){
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
find_consecutive_threshold <- function(
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
#' @param histogram_count_threshold a hard threshold to filter histogram density
#' @param eps numeric hyperparameter to finetune segmentation. See `Delon et al, 2005`
#' @param seed numeric seed
#' @param truncated_models logical, whether to fit truncated distributions
#' @param max_uniform logical, whether to find a subsegment maximizing the fit of a uniform distribution
#' @param uniform_peak_threshold numeric, indicating the minimum proportion of the subsegment
#' @param uniform_peak_stepsize integer, indicating the stepsize (relative to the histogram bins) to take in the search for the uniform subsegment
#' @param remove_low_entropy logical, indicating whether to filter out low entropy regions
#' @param min_gap_size integer, indicating the minimum gap size to be filtered
#' @param histogram_metric a subset of `jaccard`, `intersection`, `ks`, `mse`, `chisq` indicating metrics to use for fit optimization
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

  # Error checking
  stopifnot(inherits(histogram_obj, "Histogram"))
  if(!is.numeric(histogram_count_threshold) | !(histogram_count_threshold >= 0)){
    stop("histogram_count_threshold should be a numeric greater than or equal to 0")
  }
  if(!is.numeric(eps) | !(eps > 0)){
    stop("eps must be a positive numeric")
  }
  if(!is.logical(truncated_models) | length(truncated_models) != 1){
    stop("truncated_models has to be a logical of length 1")
  }
  if(!is.numeric(uniform_peak_threshold) | length(uniform_peak_threshold) != 1 ){
    stop("uniform_peak_threshold must be a numeric of length 1")
  }
  if(!(uniform_peak_threshold >= 0 & uniform_peak_threshold <= 1)) {
    stop("uniform_peak_threshold must be between 0 and 1")
  }
  if(!is_equal_integer(uniform_peak_stepsize) | !(uniform_peak_stepsize > 0) | length(uniform_peak_stepsize) != 1){
    stop("uniform_peak_stepsize must be functional as a positive integer of length 1")
  }
  if(!is.logical(remove_low_entropy) | length(remove_low_entropy) != 1){
    stop("remove_low_entropy has to be a logical of length 1")
  }
  if(!is_equal_integer(min_gap_size) | length(min_gap_size) != 1){
    stop("min_gap_size must be a numeric integer of length 1")
  }
  if(!is_equal_integer(min_peak_size) | length(min_peak_size) != 1){
    stop("min_peak_size must be a numeric integer of length 1")
  }
  if(!is.logical(max_uniform) | length(max_uniform) != 1){
    stop("max_uniform has to be a logical of length 1")
  }
  histogram_metric <- match.arg(histogram_metric, several.ok = T)
  distributions <- match.arg(distributions, several.ok = T)

  # Extracting data
  x <- histogram_obj$histogram_data

  # Change points
  chgpts <- find_stepfunction_chgpts(x)
  # Looking for regions that surpass a hard count threshold
  x_segs <- as.data.frame(find_consecutive_threshold(x, threshold = histogram_count_threshold))
  x_segs <- x_segs[x_segs$start != x_segs$end,]

  all_points <- apply(x_segs, 1, function(segs) {
    p_init <- unname(c(segs['start'], chgpts[chgpts > segs['start'] & chgpts < segs['end']], segs['end']))
    p_init <- sort(unique(p_init)) # meaningful gaps local also needs p_init to be sorted so temporarily adding this back
    p <- ftc(x, p_init, eps)

    # Max Gap
    if(remove_low_entropy) {
      mgaps <-  meaningful_gaps_local(x = x, seg_points = p, change_points = p_init, min_gap = min_gap_size)
      p <- p[(abs(p - segs['start']) > min_peak_size & abs(p - segs['end']) > min_peak_size) | p %in% segs]
      p_pairs <- remove_max_gaps(start_end_points = index_to_start_end(p), max_gaps = mgaps, remove_short_segment = min_peak_size) # remove.short.segment can also be used to filter min_peak_size, but doesn't extend to non remove low entropy cases
    } else {
      p <- p[(abs(p - segs['start']) > min_gap_size & abs(p - segs['end']) > min_peak_size) | p %in% segs]
      p_pairs <- index_to_start_end(p)
    }

    p_pairs
  })

  # Combine the results from each segment
  all_points <- do.call('rbind.data.frame', all_points)
  all_points <- all_points[(all_points$end - all_points$start + 1) > min_peak_size,, drop = FALSE] # Current solution for min_peak_size, open to alternatives
  rownames(all_points) <- NULL # use reset.rownames?

  # Fitting different models
  models <- list()
  set.seed(seed)
  for(i in seq_len(nrow(all_points))) {
    seg <- all_points[i, ]
    seg_start <- seg[['start']]
    seg_end <- seg[['end']]
    seg_len <- seg_end - seg_start + 1
    bin_data <- x[seg_start:seg_end]

    dist_optim <- fit_distributions(bin_data, metric = histogram_metric, truncated = truncated_models, distributions = distributions)
    dist_optim <- lapply(dist_optim, function(y) {
      y$seg_start <- seg_start
      y$seg_end <- seg_end
      y
    })

    # Find the maximum uniform segment
    if(max_uniform & seg_len > uniform_peak_stepsize & seg_len > ceiling(uniform_peak_threshold*seg_len)){
      unif_segment <- identify_uniform_segment(bin_data, metric = histogram_metric, threshold = uniform_peak_threshold, stepsize = uniform_peak_stepsize, max_sd_size = 0)
      # Use the maximum segment
      bin_data_subset <- bin_data[unif_segment[1, "start"]:unif_segment[1, "end"]]
      # Fit uniform distribution on maximum uniform segment
      dist_optim_subset <- fit_distributions(bin_data_subset, metric = histogram_metric, truncated = FALSE, distributions = "unif")
      # Adjust the segment starts from the shifted max uniform segment
      dist_optim_subset <- lapply(dist_optim_subset, function(y) {
        y$seg_start <- unif_segment[1, "start"] + seg_start - 1
        y$seg_end <- unif_segment[1, "end"] + seg_start - 1
        y$dist <- "unif"
        y
      })
      for(munif in names(dist_optim_subset)){dist_optim[[munif]] <- dist_optim_subset[[munif]]}
    }

    # Metric Voting
    metric_idx <- sapply(dist_optim, `[[`, "metric")
    best_models <- lapply(histogram_metric, function(met){
      metric_histograms <- dist_optim[which(metric_idx == met)]
      metric_fit <- lapply(metric_histograms, `[[`, "value")
      return( metric_histograms[which.min(metric_fit)] )
    })
    best_models <- unlist(best_models, recursive = F)
    dist_extract <- sapply(best_models, `[[`, "dist")
    dist_vote <- sort(table(dist_extract), decreasing = T)
    if(length(dist_vote) > 1 & dist_vote[1] == dist_vote[2]){
       mets <- sapply(best_models, `[[`, "metric")
       dist_majority <- dist_extract[which(mets == "jaccard")]
    } else {
      dist_majority <- names(dist_vote)[1]
    }
    best_models[['consensus']] <- dist_optim[[paste0("jaccard.", dist_majority)]]

    # Correcting for optimization via finding the minimum
    best_models <- lapply(best_models, function(mod){
      if(mod$metric %in% c("jaccard", "intersection")){ mod$value = 1 - mod$value }
      mod
    })

    best_models[['consensus']]$metric <- "consensus"
    names(best_models) <- gsub("[.].*", "", names(best_models))


    models[[i]] <- best_models
  }

  # Creating a HistogramFit object
  res <- list("models" = models, "p" = all_points,  "histogram_count_threshold" = histogram_count_threshold,
             "eps" =  eps, "seed" = seed, "truncated_models" = truncated_models, "uniform_peak_threshold" = uniform_peak_threshold,
             "uniform_peak_stepsize" = uniform_peak_stepsize, "remove_low_entropy" = remove_low_entropy, "min_gap_size" = min_gap_size,
             "min_peak_size" = min_peak_size, "max_uniform" = max_uniform, "histogram_metric" = histogram_metric, "distributions" = distributions)
  res <- c(histogram_obj, res)
  class(res) <- c("HistogramFit", class(histogram_obj))

  return(res)
}
