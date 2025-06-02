#!/bin/bash
gcc -o sequential sequential.c -lm
./sequential 512 200 1000 > sequential-512.txt