#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output-1024.txt
#SBATCH -e output/errors-1024.txt

./sequential 1024 200 1000

