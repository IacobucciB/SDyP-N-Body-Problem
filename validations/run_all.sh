#!/bin/bash

# Limpiar archivos de salida previos
> all_validations_pthreads.txt
> all_validations_mpi.txt

# Copiar y renombrar archivos de salida secuencial en todas las carpetas de Pthreads y openMPI-Pthreads
for size in 512 1024 2048 4096; do
  src_file="../scripts/Secuencial/output/output-${size}.txt"
  # Pthreads
  for core in 4-cores 8-cores; do
    dest_dir="Pthreads/$core/$size"
    mkdir -p "$dest_dir"
    if [ -f "$src_file" ]; then
      cp "$src_file" "$dest_dir/ARCHIVO-SALIDA-SECUENCIAL.txt"
    fi
  done
  # openMPI-Pthreads
  for core in 4-cores 8-cores 16-cores; do
    dest_dir="openMPI-Pthreads/$core/$size"
    mkdir -p "$dest_dir"
    if [ -f "$src_file" ]; then
      cp "$src_file" "$dest_dir/ARCHIVO-SALIDA-SECUENCIAL.txt"
    fi
  done
done

for core in 4-cores 8-cores; do
  for size in 512 1024 2048 4096; do
    # Copiar y renombrar archivo de salida a validations
    src_file="../scripts/Pthreads/$core/output/output${size}.txt"
    dest_dir="Pthreads/$core/$size"
    mkdir -p "$dest_dir"
    if [ -f "$src_file" ]; then
      cp "$src_file" "$dest_dir/ARCHIVO-SALIDA-PARALELO.txt"
    fi
  done
done

# Ejecutar todos los run_validation.sh en Pthreads
echo "\n==== Pthreads ===="
for core in 4-cores 8-cores; do
  for size in 512 1024 2048 4096; do
    script="Pthreads/$core/$size/run_validation.sh"
    if [ -f "$script" ]; then
      echo -e "\n==== Running $script ====" >> all_validations_pthreads.txt
      bash "$script" >> all_validations_pthreads.txt
    fi
  done
done

for core in 4-cores 8-cores 16-cores; do
  for size in 512 1024 2048 4096; do
    # Copiar y renombrar archivo de salida a validations
    src_file="../scripts/openMPI-Pthreads/$core/output/output${size}.txt"
    dest_dir="openMPI-Pthreads/$core/$size"
    mkdir -p "$dest_dir"
    if [ -f "$src_file" ]; then
      cp "$src_file" "$dest_dir/ARCHIVO-SALIDA-PARALELO.txt"
    fi
  done
done

# Ejecutar todos los run_validation.sh en openMPI-Pthreads
echo "\n==== openMPI-Pthreads ===="
for core in 4-cores 8-cores 16-cores; do
  for size in 512 1024 2048 4096; do
    script="openMPI-Pthreads/$core/$size/run_validation.sh"
    if [ -f "$script" ]; then
      echo -e "\n==== Running $script ====" >> all_validations_mpi.txt
      bash "$script" >> all_validations_mpi.txt
    fi
  done
done