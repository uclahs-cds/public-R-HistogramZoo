#!/bin/bash

while true
do
  for i in {1..12}
  do
    sbatch sbatch_unimodal.sh 100
  done
  sleep 1000
  echo "New loop"
done
