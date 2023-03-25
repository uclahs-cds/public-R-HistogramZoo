#!/bin/sh
#SBATCH --partition=F2

######################
echo $1
Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/unimodal_simulation_param_sample.R $1
