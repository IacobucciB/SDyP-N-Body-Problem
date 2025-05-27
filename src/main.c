#include <stdio.h>
#include "pthreads_source.h"

int main() {
    printf("Hello, world!\n");

    int rank = 4; // Example rank, can be set dynamically
    pthreads_functions(rank);

    return 0;
}
