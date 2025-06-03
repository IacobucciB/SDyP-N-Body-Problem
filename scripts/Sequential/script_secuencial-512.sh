#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output-512.txt
#SBATCH -e output/errors-512.txt

./n_body_simple_NOGL 512 200 1000

