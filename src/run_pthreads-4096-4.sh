#!/bin/bash 
gcc -o pthreads pthreads.c -pthread -lm
./pthreads 4096 200 1000 4 > pthreads-4096-4.txt

