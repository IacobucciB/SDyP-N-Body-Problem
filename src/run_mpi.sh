#!/bin/bash

mpicc -o mpi mpi_new.c -lm -pthread
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi

mpirun --oversubscribe --bind-to none -np 2 mpi 512 200 1000 2 > output_mpi.txt
if [ $? -ne 0 ]; then
    echo "Error: MPI execution failed."
    exit 1
fi

echo "Execution completed successfully."