#!/bin/bash

# Run all other run_validation.sh scripts in Pthreads subfolders
for core in 4-cores 8-cores; do
  for size in 1024 512; do
    script="Pthreads/$core/$size/run_validation.sh"
    if [ -f "$script" ]; then
      echo "\n==== Running $script ====" >> all_validations_pthreads.txt
      bash "$script" >> all_validations_pthreads.txt
    fi
  done
done

# Run all other run_validation.sh scripts in openMPI-Pthreads subfolders
for core in 4-cores 8-cores 16-cores; do
  for size in 1024 2048 4096 512; do
    script="openMPI-Pthreads/$core/$size/run_validation.sh"
    if [ -f "$script" ]; then
      echo "\n==== Running $script ====" >> all_validations_mpi.txt
      bash "$script" >> all_validations_mpi.txt
    fi
  done
done