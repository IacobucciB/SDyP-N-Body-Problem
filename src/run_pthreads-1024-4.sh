#!/bin/bash 
gcc -o pthreads pthreads.c -pthread -lm
./pthreads 1024 200 1000 4 > pthreads-1024-4.txt
