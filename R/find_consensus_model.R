
#' Methods for voting for a consensus model based on the metrics of fit_distributions
#'
#' @param models a list of models, e.g. derived from `fit_distributions`
#' @param method one of `weighted_majority_vote` and `rra` as a method of determining the best method
#' @param metrics metrics used to fit models. Metrics should be ordered in descending priority. The first metric in the vector will be used to return the `consensus` model for the distribution determined through voting.
#' @param weights required if `method` is `weighted_majority_voting`. weights of each metric to be multiplied by rankings. Weights should be in decreasing order. A higher weight results in a higher priority of the metric.
#'
#' @return a list of the best model for each metric and a `consensus` model representing the model with the consensus distribution and the lowest weighted metric.
#' 
#' @export
#'
#' @examples \dontrun{
#'  data <- observations_to_histogram(rnorm(10000, mean = 20, sd = 10))
#'  data <- data$histogram_data
#'  models <- fit_distributions(data)
#'  find_consensus_models(models)
#' 
#' }
find_consensus_model <- function(
    models,
    method = c("weighted_majority_vote", "rra"),
    metrics = c("jaccard", "intersection", "ks", "mse", "chisq"),
    weights = rev(seq(1, 1.8, 0.2))[1:length(metrics)]
){

  # Initialization
  met <- sapply(models, `[[`, "metric")
  dist <- sapply(models, `[[`, "dist")
  val <- sapply(models, `[[`, "value")
  tag <- paste0(met, ".", dist)
  
  # Error checking
  metrics <- match.arg(metrics, several.ok = T)
  method <- match.arg(method)
  stopifnot(is.numeric(weights))
  if(!all(sort(weights, decreasing = T) == weights)){
    warning("Weights should be in decreasing order.")
  }
  if(length(weights) != length(metrics) & method == "weighted_majority_voting"){
    stop("Numeric weights must be provided for all metrics.")
  }
  if(any(sapply(met, is.null))){
    stop("One or more models are missing metrics.")
  }
  if(!all(metrics %in% met)){
    stop("One or more provided metrics are not found in models.")
  }
  if(any(sapply(dist, is.null))){
    stop("One or more models are missing distributions.")
  }
  if(any(sapply(val, is.null))){
    stop("One or more models are missing fitted values.")
  }
  
  # Set-up
  model_metrics <- cbind.data.frame(met, dist, val)
  model_metrics <- reshape2::acast(
    data = model_metrics, 
    formula = dist ~ met, 
    value.var = "val"
  )
  model_metrics <- model_metrics[,metrics]
  
  # Majority voting
  if(method == "weighted_majority_vote"){
    model_metrics[] <- apply(model_metrics, 2, function(x) rank(-x))
    model_metrics[] <- model_metrics %*% diag(weights)
    score <- rowSums(model_metrics)
    if(sum(score == max(score)) > 1){
      stop("Ties exist between distributions chosen by metrics.")
    }
    model_metrics <- model_metrics[order(score, decreasing = T),]
    consensus_dist <- rownames(model_metrics)[1]
    metric_best_model <- apply(model_metrics, 2, which.max)
  }
  
  # Rank aggregation
  if(method == "rra"){
    model_metrics[] <- apply(model_metrics, 2, function(x) rank(x)/length(x))
    rra <- RobustRankAggreg::aggregateRanks(
      rmat = model_metrics, 
      N = nrow(model_metrics)
    )
    if(sum(rra$Score == min(rra$Score)) > 1){
      stop("Ties exist between distributions chosen by metrics.")
    }
    consensus_dist <- rra[1, "Name"]
    metric_best_model <- apply(model_metrics, 2, which.min)
  }
  
  # Extracting the best models
  best_model_tags <- paste0(metrics, ".", rownames(model_metrics)[metric_best_model])
  best_models <- models[tag %in% best_model_tags]
  names(best_models) <- sapply(best_models, `[[`, "metric")
  best_models[['consensus']] <- models[[which(tag == paste0(metrics[1], ".", consensus_dist))]]
  
  return(
    best_models
  )
}