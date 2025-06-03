#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output4096.txt
#SBATCH -e output/errors4096.txt
./pthreads 4096 200 1000 4

