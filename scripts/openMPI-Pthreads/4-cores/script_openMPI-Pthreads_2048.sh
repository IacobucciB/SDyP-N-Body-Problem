#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o output/output2048.txt
#SBATCH -e output/errors2048.txt
mpirun --bind-to none -np 2 ./mpi 2048 200 1000 2

