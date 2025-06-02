#!/bin/bash
gcc -o sequential sequential.c -lm
./sequential 4096 200 1000 > sequential-4096.txt

