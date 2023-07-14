#!/bin/sh
#SBATCH --partition=F2
#SBATCH --mem-per-cpu=1gb
#SBATCH --ntasks=1
#SBATCH -n 1
##SBATCH --exclusive
set -e
######################
echo $1

Rscript -e 'remotes::install_local("/hot/users/stefaneng/public-R-HistogramZoo", dependencies=FALSE)'
# R CMD install /hot/users/stefaneng/public-R-HistogramZoo

while true
do
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/unimodal_simulation_param_sample.R $1
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/unimodal_simulation_param_sample.R $1 MLE
done
