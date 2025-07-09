#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o output/output512.txt
#SBATCH -e output/errors512.txt
mpirun --bind-to none -np 2 ./mpi 512 200 1000 4

