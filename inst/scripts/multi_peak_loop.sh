#!/bin/sh
#SBATCH --partition=F2
#SBATCH --mem-per-cpu=1gb
#SBATCH --ntasks=1
#SBATCH -n 1

set -e
######################
echo $1

# R CMD install /hot/users/stefaneng/public-R-HistogramZoo
Rscript -e 'remotes::install_local("/hot/users/stefaneng/public-R-HistogramZoo", dependencies=FALSE)'

while true
do
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/multimodal_simulation_param_sample.R $1
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/multimodal_simulation_param_sample.R $1 MLE
done
