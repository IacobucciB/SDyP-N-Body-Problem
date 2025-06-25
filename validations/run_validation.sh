#!/bin/bash

# Limpiar archivos de salida previos
> all_validations_pthreads.txt
> all_validations_mpi.txt

# Ejecutar todos los run_validation.sh en Pthreads
for core in 4-cores 8-cores; do
  for size in 512 1024 2048 4096; do
    script="Pthreads/$core/$size/run_validation.sh"
    if [ -f "$script" ]; then
      echo -e "\n==== Running $script ====" >> all_validations_pthreads.txt
      bash "$script" >> all_validations_pthreads.txt
    fi
  done
done

# Ejecutar todos los run_validation.sh en openMPI-Pthreads
for core in 4-cores 8-cores 16-cores; do
  for size in 512 1024 2048 4096; do
    script="openMPI-Pthreads/$core/$size/run_validation.sh"
    if [ -f "$script" ]; then
      echo -e "\n==== Running $script ====" >> all_validations_mpi.txt
      bash "$script" >> all_validations_mpi.txt
    fi
  done
done