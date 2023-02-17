library(HistogramZoo)

output_path <- '/hot/users/stefaneng/public-R-HistogramZoo/results'
set.seed(314)

# Parameters that we want to vary:
# Sample size
# Normal dist (mean, sd)
# Skew (gamma case)
# Number of distributions
# Time

common_params <- list(
  N = c(1e4, 1e5, 1e6),
  seed = sample(1:1e5, size = 20)
)

norm_params <- do.call('expand.grid', c(
  common_params,
  list(
      param1 = 50,
      param2 = seq(5, 100, length.out = 15)
    )
  ))

norm_args <- list(
    eps = 1,
    remove_low_entropy = TRUE,
    truncated_models = FALSE
  )

cat('Running unimodal norm simulation:\n')
norm_results <- do.call(
  'mapply',
  c(list(FUN = HistogramZoo:::general_sim, SIMPLIFY = FALSE, quiet = F, rfunc = 'rnorm'),
    norm_args,
    norm_params)
  )

saveRDS(norm_results, file = file.path(output_path, 'norm_sim_unimodal.rds'))

## Gamma sim
gamma_params <- do.call('expand.grid', c(
    common_params,
    list(
        param1 = seq(0.4, 25, length.out = 25),
        param2 = 1
      )
    )
  )

gamma_args <- list(
    eps = 1,
    remove_low_entropy = FALSE,
    truncated_models = FALSE
  )

cat('Running unimodal gamma simulation:\n')
gamma_results <- do.call(
  'mapply',
  c(list(FUN = HistogramZoo:::general_sim, SIMPLIFY = FALSE, quiet = F, rfunc = 'rgamma'),
    gamma_args,
    gamma_params)
  )

saveRDS(gamma_results, file = file.path(output_path, 'gamma_sim_unimodal.rds'))

### Unif
unif_params <- do.call('expand.grid', c(
    common_params,
    list(
        param1 = 0, # min
        param2 = 10 # max
      )
    )
  )

unif_args <- list(
    eps = 1,
    remove_low_entropy = FALSE,
    truncated_models = FALSE
  )

cat('Running unimodal unif simulation:\n')
unif_results <- do.call(
  'mapply',
  c(list(FUN = HistogramZoo:::general_sim, SIMPLIFY = FALSE, quiet = F, rfunc = 'runif'),
    unif_args,
    unif_params)
  )

saveRDS(norm_results, file = file.path(output_path, 'unif_sim_unimodal.rds'))

cat('Running unimodal log-normal simulation:\n')

lognorm_params <- do.call('expand.grid', c(
    common_params,
    list(
        param1 = 0,
        param2 = seq(0.1, 3, length.out = 15)
      )
    )
  )

lognorm_args <- list(
    eps = 1,
    remove_low_entropy = FALSE,
    truncated_models = FALSE
  )

cat('Running log-normal simulation:\n')
unif_results <- do.call(
  'mapply',
  c(list(FUN = HistogramZoo:::general_sim, SIMPLIFY = FALSE, quiet = F, rfunc = 'rlnorm'),
    lognorm_args,
    lognorm_params)
  )

saveRDS(norm_results, file = file.path(output_path, 'lognormal_sim_unimodal.rds'))
