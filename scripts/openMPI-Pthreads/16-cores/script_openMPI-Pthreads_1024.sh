#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o output/output1024.txt
#SBATCH -e output/errors1024.txt
mpirun --bind-to none mpi 1024 200 1000 8

