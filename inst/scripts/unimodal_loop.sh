#!/bin/bash

while true
do
  for i in {1..20}
  do
    sbatch sbatch_unimodal.sh 300
  done 
  sleep 900
  echo "New loop"
done
