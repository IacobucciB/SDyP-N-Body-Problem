#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o output/output.txt
#SBATCH -e output/errors.txt

./n_body_simple_NOGL 512 200 1000
