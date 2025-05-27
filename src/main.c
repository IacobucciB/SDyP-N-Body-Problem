// Compile:
//      gcc -o main src/main.c src/pthreads_source.c -lpthread
// Run:
//      ./main

#include <stdio.h>
#include "pthreads_source.h"

int main(int argc, char *argv[]){
    printf("Hello, world!\n");

    if (argc < 4) {
        printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos>\n", argv[0]);
        return -1;
    }

    int rank = 4; // Example rank, can be set dynamically
    pthreads_functions(rank);

    return 0;
}
