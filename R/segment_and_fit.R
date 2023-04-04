#' Find changepoints in a vector with uniform stretches of values
#' @param x A numeric vector
find_change_points <- function(x){
  change_points <- which(diff(x) != 0)
  change_points_plus <- change_points+1
  keep <- (x[change_points] < x[change_points_plus])
  return( c(change_points[keep], change_points_plus[!keep]) )
}


#' Returns the indices for consecutive elements of a vector that are greater than a specified threshold
#'
#' @param x numeric vector
#' @param threshold numeric threshold
#'
#' @return list of coordinates with `start` and `end` coordinates
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

#' Segmentation of histograms and distribution fitting
#'
#' @param histogram_obj a Histogram or HistogramList object
#' @param optima_threshold threshold for local optima, i.e. a point can only be considered a local optima if it differs from its neighbour optima by greater than the permitted threshold, default 0
#' @param optima_flat_endpoints in regions of flat density, whether to return the endpoints or the midpoints
#' @param histogram_count_threshold a hard threshold to filter histogram density
#' @param eps numeric (epsilon) hyperparameter to finetune segmentation. See `Delon et al, 2005`
#' @param remove_low_entropy logical, indicating whether to filter out low entropy regions
#' @param min_gap_size integer, indicating the minimum gap size to be filtered
#' @param min_segment_size integer, indication the minimum segment size, default 3
#' @param seed numeric seed
#' @param max_uniform logical, whether to find a subsegment maximizing the fit of a uniform distribution
#' @param uniform_threshold numeric, indicating the minimum proportion of the subsegment
#' @param uniform_stepsize integer, indicating the stepsize (relative to the histogram bins) to take in the search for the uniform subsegment
#' @param uniform_max_sd numeric, the number of standard deviations of the computed metric distribution away from the optimal uniform which has maximum length
#' @param truncated_models logical, whether to fit truncated distributions
#' @param metric a subset of `jaccard`, `intersection`, `ks`, `mse`, `chisq` indicating metrics to use for fit optimization. Metrics should be ordered in descending priority. The first metric in the vector will be used to return the `consensus` model for the distribution determined through voting.
#' @param distributions a subset of `norm`, `gamma`, and `unif` indicating distributions to fit.
#' @param consensus_method one of `weighted_majority_vote` and `rra` as a method of determining the best method
#' @param metric_weights required if `method` is `weighted_majority_voting`. weights of each metric to be multiplied by rankings. Weights should be in decreasing order. A higher weight results in a higher priority of the metric.
#' @param distribution_prioritization if `method` is `weighted_majority_voting`, a list of ranked distributions, to break ties
#'
#' @return a HistogramFit object representing the Histogram and results of the fit
#' @export
#' @examples \dontrun{
#' x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
#' res = segment_and_fit(x)
#' }
segment_and_fit <- function(
    histogram_obj,
    optima_threshold = 0,
    optima_flat_endpoints = T,
    histogram_count_threshold = 0,
    eps = 1,
    remove_low_entropy = T,
    min_gap_size = 2,
    min_segment_size = 3,
    seed = NULL,
    max_uniform = T,
    uniform_threshold = 0.75,
    uniform_stepsize = 5,
    uniform_max_sd = 0,
    truncated_models = FALSE,
    metric = c("mle", "jaccard", "intersection", "ks", "mse", "chisq"),
    distributions = c("norm", "unif", "gamma", "gamma_flip"),
    consensus_method = c("weighted_majority_vote", "rra"),
    metric_weights = rev(seq(1, 1.8, 0.2)),
    distribution_prioritization = distributions
) {

  # Error checking
  stopifnot(inherits(histogram_obj, "Histogram"))
  if(!is.numeric(histogram_count_threshold) | !(histogram_count_threshold >= 0)){
    stop("histogram_count_threshold should be a numeric greater than or equal to 0")
  }
  if(!is_equal_integer(optima_threshold) | !(optima_threshold >= 0) | length(optima_threshold) != 1){
    stop("optima_threshold must be functional as a integer greater than or equal to 0 and of length 1")
  }
  if(!is.logical(optima_flat_endpoints) | length(optima_flat_endpoints) != 1){
    stop("optima_flat_endpoints has to be a logical of length 1")
  }
  if(!is.numeric(eps) | !(eps > 0)){
    stop("eps must be a positive numeric")
  }
  if(!is.logical(truncated_models) | length(truncated_models) != 1){
    stop("truncated_models has to be a logical of length 1")
  }
  if(!is.numeric(uniform_threshold) | length(uniform_threshold) != 1 ){
    stop("uniform_threshold must be a numeric of length 1")
  }
  if(!(uniform_threshold >= 0 & uniform_threshold <= 1)) {
    stop("uniform_threshold must be between 0 and 1")
  }
  if(!is_equal_integer(uniform_stepsize) | !(uniform_stepsize > 0) | length(uniform_stepsize) != 1){
    stop("uniform_stepsize must be functional as a positive integer of length 1")
  }
  if(!is.numeric(uniform_max_sd) | !(uniform_max_sd >= 0) | length(uniform_max_sd) != 1){
    stop("uniform_max_sd must be a positive or zero numeric of length 1")
  }
  if(!is.logical(remove_low_entropy) | length(remove_low_entropy) != 1){
    stop("remove_low_entropy has to be a logical of length 1")
  }
  if(!is_equal_integer(min_gap_size) | length(min_gap_size) != 1){
    stop("min_gap_size must be a numeric integer of length 1")
  }
  if(!is_equal_integer(min_segment_size) | length(min_segment_size) != 1){
    stop("min_segment_size must be a numeric integer of length 1")
  }
  if(!is.logical(max_uniform) | length(max_uniform) != 1){
    stop("max_uniform has to be a logical of length 1")
  }

  metric <- match.arg(metric, several.ok = T)
  consensus_method <- match.arg(consensus_method)
  distributions <- match.arg(distributions, several.ok = T)

  # Pre-error checking for weighted_majority_vote
  if(consensus_method == "weighted_majority_vote"){
    if(!all(sort(distributions) == sort(distribution_prioritization))){
      stop("distribution_prioritization needs be an ordering of distributions")
    }
    if(!is.numeric(metric_weights) | length(metric_weights) != length(metric)){
      stop("Numeric weights must be provided for all metrics.")
    }
    if(!all(sort(metric_weights, decreasing = T) == metric_weights)){
      warning("Weights should be in decreasing order.")
    }
  }

  # Extracting data
  x <- histogram_obj$histogram_data

  # Finding local optima
  optima <- find_local_optima(x, threshold = optima_threshold, flat_endpoints = optima_flat_endpoints)
  optima <- sort(c(optima$min_ind, optima$max_ind))

  # Finding change points for remove_low_entropy
  if(remove_low_entropy){
    changepoints <- find_change_points(x)
  }

  # Looking for regions that surpass a hard count threshold
  x_segs <- as.data.frame(find_consecutive_threshold(x, threshold = histogram_count_threshold))
  x_segs <- x_segs[x_segs$start != x_segs$end,]

  if(nrow(x_segs) == 0){
    stop("no segments of greater than length 1 remain after filtering for histogram_count_threshold")
  }

  # Identifying endpoints of each segment
  all_points <- apply(x_segs, 1, function(segs) {
    p_init <- unname(c(segs['start'], optima[optima > segs['start'] & optima < segs['end']], segs['end']))
    p <- ftc(x, p_init, eps)

    # Max gap
    if(remove_low_entropy) {
      changepoints_subset <- unname(c(segs['start'], changepoints[changepoints > segs['start'] & changepoints < segs['end']], segs['end']))
      mgaps <-  meaningful_gaps_local(x = x, seg_points = p, change_points = changepoints_subset, min_gap = min_gap_size)
      p <- p[((p - segs['start'] + 1) >= min_segment_size & (segs['end'] - p + 1) >= min_segment_size) | p %in% segs]
      p_pairs <- remove_max_gaps(start_end_points = index_to_start_end(p), max_gaps = mgaps, remove_short_segment = min_segment_size)
    } else {
      p <- p[((p - segs['start'] + 1) >= min_segment_size & (segs['end'] - p + 1) >= min_segment_size) | p %in% segs]
      p_pairs <- index_to_start_end(p)
    }

    p_pairs
  })

  # Combine the results from each segment
  all_points <- do.call('rbind.data.frame', all_points)
  all_points <- all_points[(all_points$end - all_points$start + 1) >= min_segment_size,, drop = FALSE]
  rownames(all_points) <- NULL

  if(nrow(all_points) == 0){
    stop("no segments of greater than min_segment_size remain after segmentation")
  }

  # Fitting different models
  set.seed(seed)
  models <- apply(all_points, 1, function(seg){
    seg_start <- seg[['start']]
    seg_end <- seg[['end']]
    seg_len <- seg_end - seg_start + 1
    bin_data <- x[seg_start:seg_end]
    # sub_hist <- histogram_obj[seg_start:seg_end]

    # Find the maximum uniform segment
    dist_optim <- unif_segment <- list()
    if("unif" %in% distributions &&
       max_uniform &&
       seg_len > uniform_stepsize &&
       seg_len > ceiling(uniform_threshold*seg_len) &&
       ! 'mle' %in% metric
    ){
      unif_segment <- lapply(metric, function(met) {
        res <- identify_uniform_segment(
          x = bin_data,
          metric = met,
          threshold = uniform_threshold,
          stepsize = uniform_stepsize,
          max_sd_size = uniform_max_sd
        )
        res[['seg_start']] <- res[['seg_start']] + seg_start - 1
        res[['seg_end']] <- res[['seg_end']] + seg_start - 1
        return(res)
      })
      # Removing unif from the vector of distributions
      distributions <- setdiff(distributions, "unif")
    }

    if (max_uniform && 'mle' %in% metric) {
      warning("Cannot fit max_uniform with 'mle' in metrics. Setting max_uniform = FALSE.")
    }

    if( length(distributions) > 0 ){
      dist_optim <- fit_distributions(
        x = bin_data,
        metric = metric,
        truncated = truncated_models,
        distributions = distributions
      )
      dist_optim <- lapply(dist_optim, function(obj){
        obj[['seg_start']] <- seg_start
        obj[['seg_end']] <- seg_end
        return(obj)
      })
    }
    dist_optim <- c(dist_optim, unif_segment)

    # Correcting for optimization via finding the minimum
    best_models <- find_consensus_model(
      models = dist_optim,
      method = consensus_method,
      metric = metric,
      distribution_prioritization = distribution_prioritization,
      weights = metric_weights
    )

    return(best_models)
  })

  # Creating a HistogramFit object
  res <- list("models" = models, "p" = all_points,
              "optima_threshold" = optima_threshold, "optima_flat_endpoints" = optima_flat_endpoints,
              "histogram_count_threshold" = histogram_count_threshold,
              "eps" =  eps,
              "remove_low_entropy" = remove_low_entropy, "min_gap_size" = min_gap_size, "min_segment_size" = min_segment_size,
              "seed" = seed,
              "max_uniform" = max_uniform, "uniform_threshold" = uniform_threshold, "uniform_stepsize" = uniform_stepsize, "uniform_max_sd" = uniform_max_sd,
              "truncated_models" = truncated_models, "metric" = metric, "distributions" = distributions,
              "consensus_method" = consensus_method, "distribution_prioritization" = distribution_prioritization, "metric_weights" = metric_weights)
  res <- c(histogram_obj, res)
  class(res) <- c("HistogramFit", class(histogram_obj))

  return(res)
}
