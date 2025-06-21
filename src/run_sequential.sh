#!/bin/bash

gcc -o sequential sequential.c -lm
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi

./sequential 512 200 1000 > output_sequential.txt
if [ $? -ne 0 ]; then
    echo "Error: Execution failed."
    exit 1
fi

echo "Execution completed successfully."