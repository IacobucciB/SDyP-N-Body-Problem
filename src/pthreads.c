#include <pthread.h>
#include <stdio.h>

int rank;
int T;

void *pfunction(void *arg) {
    int myRank = *((int *)arg);
    printf("Hello from thread %d of %d\n", myRank, T);
    pthread_exit(NULL);
}

void pthreads_functions(int myRank) {
    rank = myRank;
    T = 4; // Example number of threads, can be set dynamically
    pthread_t threads[T];
    int threads_ids[T];
    int id;

    for (id = 0; id < T; id++) {
        threads_ids[id] = id;
        if (pthread_create(&threads[id], NULL, pfunction, (void *)&threads_ids[id]) != 0) {
            perror("Failed to create thread");
            return;
        }
    }

    for (id = 0; id < T; id++) {
        if (pthread_join(threads[id], NULL) != 0) {
            perror("Failed to join thread");
            return;
        }
    }
    

}

int main(int argc, char const *argv[])
{
    int myRank = 0; // Example rank, can be set dynamically
    pthreads_functions(myRank);
    return 0;
}
