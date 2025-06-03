#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output2048.txt
#SBATCH -e output/errors2048.txt
./pthreads 2048 200 1000 8
