covariate.col <- list(
  max.uniform = 'firebrick3',
  remove.low.entropy = 'dodgerblue2',
  N.decile = 'slateblue4',
  noise.decile = 'darkorange1',
  eps.decile = 'seagreen3',
  param.decile = 'yellowgreen',
  jaccard.decile = 'violetred3'
  )

acc.colour.scheme <- c('white', 'red');

segment_prob <- function(distribution = c('norm', 'gamma', 'unif'), params, a, b) {
  stopifnot(b >= a)
  distribution <- match.arg(distribution)
  # if (distribution == 'gamma' && a < 0 )
  cdf <- get(paste0('p', distribution))
  cdf.params.func <- function(q) {
    args <- params
    args$q <- q
    do.call(cdf, args)
  }
  cdf.params.func(b) - cdf.params.func(a)
}

common.sim.legend <- function(
    include.legends = c('params', 'distributions', 'quantiles'),
    params.to.include = c(
      'max_uniform', 'remove_low_entropy'
      ),
    cont.params.to.include = c(
      'N_decile', 'noise_decile',
      'eps_decile', 'param_decile',
      'jaccard_decile'
    )
  ) {
  include.legends <- match.arg(include.legends, several.ok = TRUE);
  params.to.include <- match.arg(params.to.include, several.ok = TRUE);
  if (! is.null(cont.params.to.include)) {
    cont.params.to.include <- match.arg(cont.params.to.include, several.ok = TRUE);
    }

  legend.list <- list()
  if ('params' %in% include.legends) {
    params.select <- match(params.to.include, c('max_uniform', 'remove_low_entropy'));

    legend.list <- c(
      legend.list,
      list(
        legend = list(
            colours = c(covariate.col$max.uniform, covariate.col$remove.low.entropy)[params.select],
            title = 'Parameters',
            labels = c('Max uniform', 'Remove low entropy')[params.select],
            size = 3,
            title.cex = 1,
            label.cex = 1,
            border = 'black'
            )
        )
      )
    }
  if ('distributions' %in% include.legends) {
    legend.list <- c(
      legend.list,
      list(
        legend = list(
            colours = distribution_colours[c('norm', 'gamma', 'unif')],
            title = 'Distribution',
            labels = c('Normal','Gamma', 'Uniform'),
            size = 3,
            title.cex = 1,
            label.cex = 1,
            border = 'black'
            )
        )
      )
    }

  if (! is.null(cont.params.to.include)) {
    cont.params <- gsub('_', '.', cont.params.to.include)

    cont.params.legend <- lapply(cont.params, function(p) {
    cont.param.name <- sub('[.]decile', '', p);

      list(
        colours = c('white', covariate.col[[p]]),
        labels = if (cont.param.name == 'jaccard') c('0', '1') else as.character(unimodal.params[[cont.param.name]]),
        at = c(0, 100),
        title = cont.param.name,
        angle = -90,
        width = 6,
        continuous = TRUE
        )
      })

    names(cont.params.legend) <- rep('legend', length(cont.params.legend))

    legend.list <- c(legend.list, cont.params.legend)
    }

  BoutrosLab.plotting.general::legend.grob(
    legends = legend.list,
    label.cex = 1.5,
    title.cex = 1.5,
    size = 3,
    title.just = 'left',
    layout = c(1,length(legend.list))
    )
  }



assign.cov.factor.col <- function(x, ramp.colors) {
  x.levels <- levels(x);
  k <- length(x.levels);

  colormap <- setNames(
    colorRampPalette(ramp.colors)(k),
    x.levels
    );

  colormap[x];
}

sim.plot.heatmap.cov <- function(cov.data) {
  if ('max_uniform' %in% colnames(cov.data)) {
    cov.data$max_uniform <- ifelse(cov.data$max_uniform, covariate.col$max.uniform, 'white');
  }
  if ('remove_low_entropy' %in% colnames(cov.data)) {
    cov.data$remove_low_entropy <- ifelse(cov.data$remove_low_entropy, covariate.col$remove.low.entropy, 'white');
  }
  if ('actual_dist' %in% colnames(cov.data)) {
    cov.data$actual_dist <- unname(distribution_colours[as.character(cov.data$actual_dist)]);
  }
  if ('dist' %in% colnames(cov.data)) {
    # TODO: Combine
    cov.data$dist <- unname(distribution_colours[as.character(cov.data$dist)]);
  }
  if ('jaccard_decile' %in% colnames(cov.data)) {
    cov.data$jaccard_decile <- assign.cov.factor.col(cov.data$jaccard_decile, c('white', covariate.col$jaccard.decile));
  }
  if ('N_decile' %in% colnames(cov.data)) {
    cov.data$N_decile <- assign.cov.factor.col(cov.data$N_decile, c('white', covariate.col$N.decile));
  }
  if ('noise_decile' %in% colnames(cov.data)) {
    cov.data$noise_decile <- assign.cov.factor.col(cov.data$noise_decile, c('white', covariate.col$noise.decile));
  }
  if ('eps_decile' %in% colnames(cov.data)) {
    cov.data$eps_decile <- assign.cov.factor.col(cov.data$eps_decile, c('white', covariate.col$eps.decile));
  }
  if ('param_decile' %in% colnames(cov.data)) {
    cov.data$param_decile <- assign.cov.factor.col(cov.data$param_decile, c('white', covariate.col$param.decile));
  }

  if (ncol(cov.data) == 1) cov.data <- cbind(cov.data[, 1], cov.data[, 1]);

  create.heatmap(
    cov.data,
    same.as.matrix = TRUE,
    input.colours = TRUE,
    clustering.method = 'none',
    print.colour.key = FALSE,
    xaxis.tck = 0,
    yaxis.tck = 0
    )
}
