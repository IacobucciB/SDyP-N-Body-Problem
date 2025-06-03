#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output-2048.txt
#SBATCH -e output/errors-2048.txt

./n_body_simple_NOGL 2048 200 1000

