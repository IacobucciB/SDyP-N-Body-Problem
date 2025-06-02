#!/bin/bash
gcc -o sequential sequential.c -lm
./sequential 2048 200 1000 > sequential-2048.txt

