#!/bin/sh
#SBATCH --partition=F2
### #SBATCH --nodelist=F16-2
#SBATCH -n 1
#SBATCH --mem-per-cpu=1gb
#SBATCH --ntasks=1

######################
echo $1
Rscript /hot/users/stefaneng/public-R-HistogramZoo/inst/unimodal_simulation_param_sample.R $1
