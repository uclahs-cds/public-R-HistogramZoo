
#' Methods for voting for a consensus model based on the metrics of fit_distributions
#'
#' @param models A list of models, e.g. derived from `fit_distributions`
#' @param method One of `weighted_majority_voting` and `rra` as a method of determining the best method
#' @param metrics Metrics used to fit models
#' @param inverse_weights Required if `method` is `weighted_majority_voting`. 
#' Weights of each metric to be multiplied by rankings. A lower weight results in a higher priority of lower rankings.
#'
#' @return A list of the best model for each metric and a `consensus` model representing the model with the consensus distribution and the lowest weighted metric.
#' @export
#'
#' @examples
find_consensus_model <- function(
    models,
    method = c("weighted_majority_voting", "rra"),
    metrics = c("jaccard", "intersection", "ks", "mse", "chisq"),
    inverse_weights = seq(1, 1.8, 0.2)
){

  # Error checking
  # TODO: check that each model has metrics, dist, and value
  metrics <- match.arg(metrics, several.ok = T)
  method <- match.arg(method)
  stopifnot(is.numeric(weights))
  if(length(weights) != length(metrics) & method == "weighted_majority_voting"){
    stop("Numeric weights must be provided for all metrics.")
  }
    
  # Base case - there is only one metric
  if(length(metrics) == 1){
    # pass
  }
  
  # Set-up
  dist_optim_metrics <- do.call(
    "rbind.data.frame",
    lapply(models, `[`, c("metric", "dist", "value"))
  )
  dist_optim_metrics <- reshape2::acast(
    data = dist_optim_metrics, 
    formula = dist ~ metric, 
    value.var = "value"
  )
  dist_optim_metrics[] <- apply(
    dist_optim_metrics, 2, function(x) order(x, decreasing = T)
  )
  dist_optim_metrics <- dist_optim_metrics[,metrics]
  
  # Majority voting
  if(method == "weighted_majority_voting"){
    dist_optim_metrics <- dist_optim_metrics %*% diag(weights)
    dist_optim_metrics <- dist_optim_metrics[order(rowSums(dist_optim_metrics)),]
    consensus_dist <- rownames(dist_optim_metrics)[1]
  }
  
  # Rank aggregation
  if(method == "rra"){
    dist_optim_metrics <- dist_optim_metrics/nrow(dist_optim_metrics)
    rra <- RobustRankAggreg::aggregateRanks(
      rmat = dist_optim_metrics, 
      N = nrow(dist_optim_metrics)
    )
    consensus_dist <- rra[1, "Name"]
  }
  
  # Extracting the best models
  metric_best_model <- apply(dist_optim_metrics, 2, which.min)
  best_model_names <- paste0(rownames(dist_optim_metrics)[metric_best_model], ".", metrics)
  best_models <- models[best_model_names] # What if models aren't named
  best_models[['consensus']] <- models[[paste0(metrics[1], consensus_dist)]]
  
  return(
    best_models
  )
}
