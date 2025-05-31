#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output512.txt
#SBATCH -e output/errors512.txt
./pthreads 512 200 1000 4
