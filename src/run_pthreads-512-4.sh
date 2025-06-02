#!/bin/bash
gcc -o pthreads pthreads.c -pthread -lm
./pthreads 512 200 1000 4 > pthreads-512-4.txt