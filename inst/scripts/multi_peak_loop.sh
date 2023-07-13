#!/bin/sh
#SBATCH --partition=F2
##SBATCH --exclusive

######################
echo $1

R CMD install /hot/users/stefaneng/public-R-HistogramZoo

while true
do
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/multimodal_simulation_param_sample.R $1
  Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/multimodal_simulation_param_sample.R $1 MLE
done
