#!/bin/bash

while true
do
  for i in {1..15}
  do
    sbatch sbatch_unimodal.sh 100
  done
  sleep 1000
  echo "New loop"
done
