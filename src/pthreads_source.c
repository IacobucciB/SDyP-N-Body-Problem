#include <pthread.h>
#include <stdio.h>

int rank;
int T;

void *function(void *arg) {
    int myRank = *((int *)arg);
    printf("Hello from thread %d of %d\n", myRank, T);
    return NULL;
}

void pthreads_functions(int myRank) {
    rank = myRank;
    T = 4; // Example number of threads, can be set dynamically

    pthread_t threads[T];
    int thread_args[T];

    for (int i = 0; i < T; i++) {
        thread_args[i] = i;
        pthread_create(&threads[i], NULL, function, (void *)&thread_args[i]);
    }

    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }
}