#!/bin/bash
set -e

# Rscript -e 'remotes::install_local("/hot/users/stefaneng/public-R-HistogramZoo", dependencies=FALSE, force = TRUE)'
R CMD install /hot/users/stefaneng/public-R-HistogramZoo

TODAY_DATE=$(date +'%Y-%m-%d')

# Merge the simulation data and save
Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/merge_sim_data.R
# Load the merged simulation data
Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/unimodal_analyze_sim.R ${TODAY_DATE}
