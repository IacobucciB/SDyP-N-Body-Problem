#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o output/output4096.txt
#SBATCH -e output/errors4096.txt
mpirun --bind-to none -np 2 ./mpi 4096 200 1000 4

