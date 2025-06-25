#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output-4096.txt
#SBATCH -e output/errors-4096.txt

./sequential 4096 200 1000

