#!/bin/sh
#SBATCH --partition=F2
#SBATCH --mem-per-cpu=1gb
#SBATCH --ntasks=1
#SBATCH -n 1

set -e
######################
N=${1:-50}
echo "Multi-peak sim - Number of sims per run: ${N}"

# R CMD install /hot/users/stefaneng/public-R-HistogramZoo
Rscript -e 'remotes::install_local("/hot/users/stefaneng/public-R-HistogramZoo", dependencies=FALSE)'

while true
do
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/multimodal_simulation_param_sample.R $N
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/multimodal_simulation_param_sample.R $N MLE
done
