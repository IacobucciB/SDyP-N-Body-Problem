#!/bin/bash
gcc -o sequential sequential.c -lm
./sequential 1024 200 1000 > sequential-1024.txt
