#!/bin/bash 
gcc -o pthreads pthreads.c -pthread -lm
./pthreads 2048 200 1000 4 > pthreads-2048-4.txt
