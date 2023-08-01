#!/usr/bin/env Rscript

library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');
plots.folder <- file.path(base.path, 'plots');

resolution <- 1600;

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1)

sim.version <- c('v3');
sim.merge.date <- if (length(args) == 0) Sys.Date() else args[1];

sim.version.regex <- paste0('[', paste0(sim.version, collapse = '|'), ']')
sim.version <- paste0(sim.version, collapse = '-')

overlap.files <- list.files(
  path = merged.folder,
  pattern = paste0(sim.merge.date, '_.*', sim.version.regex, '.*overlap'),
  full.names = TRUE
  )

overlap <- rbindlist(
  lapply(overlap.files, data.table::fread, sep = '\t'),
  fill = TRUE
  )

overlap.jaccard <- overlap[
  ,
  c(
    HistogramZoo:::compute_overall_jaccard(.SD),
    N = sum(.SD[fit_peak_num == 1, 'N']),
    peak_count = max(peak_num),
    fit_peak_count = max(fit_peak_num)
  ),
  by = .(seed, eps, noise, max_uniform, remove_low_entropy, mle)
  ]

print(head(overlap.jaccard))
# TODO: Heatmap of jaccard like unimodal but with peak_count
# Major differences are:
#  - No distributions
#  - N is different for each peak count = 2, 3, 4. Max is 500 * peak_count
