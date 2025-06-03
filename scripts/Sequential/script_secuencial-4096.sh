#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output-4096.txt
#SBATCH -e output/errors-4096.txt

./n_body_simple_NOGL 4096 200 1000

