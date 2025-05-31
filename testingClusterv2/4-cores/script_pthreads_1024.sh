#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output1024.txt
#SBATCH -e output/errors1024.txt
./pthreads 1024 200 1000 4

