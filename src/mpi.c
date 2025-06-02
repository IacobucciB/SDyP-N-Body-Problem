#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <mpi.h>

/* Para tiempo de ejecucion */
double dwalltime();
double tIni, tFin, tTotal;

/* Pthreads */
void *pfunction(void *arg);

/* Argumentos */
int N;        // Número de cuerpos
int dt;       // Intervalo de tiempo, longitud de un paso
int pasos;    // Número de pasos a simular
int T_MPI;    // Número total de procesos MPI
int T_PTHREADS; // Número de threads por proceso
int idW_MPI; // Identificador del proceso MPI
/* Global array */
int array[512];

pthread_barrier_t barrier; // Barrera global

/* MAIN */
int main(int argc, char *argv[])
{
    // MPI
    int mpi_error = MPI_Init(&argc, &argv); // Check MPI initialization
    if (mpi_error != MPI_SUCCESS) {
        fprintf(stderr, "Error initializing MPI.\n");
        return -1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &idW_MPI); // Obtiene el identificador de cada proceso (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &T_MPI);   // Obtiene el número total de procesos (size)

    // ARGUMENTOS
    if (argc < 5)
    {
        fprintf(stderr, "Ejecutar: %s <nro. de cuerpos> <DT> <pasos> <threads>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    N = atoi(argv[1]);
    dt = atoi(argv[2]);
    pasos = atoi(argv[3]);
    T_PTHREADS = atoi(argv[4]);

    if (N <= 0 || dt <= 0 || pasos <= 0 || T_PTHREADS <= 0) {
        fprintf(stderr, "Error: All arguments must be positive integers.\n");
        MPI_Finalize();
        return -1;
    }

    printf("N: %d, DT: %d, PASOS: %d, T_PTHREADS: %d\n", N, dt, pasos, T_PTHREADS);

    int slice_MPI = N / T_MPI; // Número de cuerpos por proceso
    int ini_MPI = idW_MPI * slice_MPI; // Índice inicial para este proceso
    int lim_MPI = ini_MPI + slice_MPI; // Índice final para este proceso

    // Inicializar el arreglo global
    for (int i = 0; i < 512; i++) {
        array[i] = i;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Cada proceso imprime su porción del arreglo
    printf("Proceso %d: Porción del arreglo [%d, %d):\n", idW_MPI, ini_MPI, lim_MPI);
    for (int i = ini_MPI; i < lim_MPI && i < 512; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");

    // Inicializar la barrera para T_PTHREADS hilos
    pthread_barrier_init(&barrier, NULL, T_PTHREADS);

    pthread_t threads[T_PTHREADS];
    int thread_ids[T_PTHREADS];
    for (int i = 0; i < T_PTHREADS; i++) {
        thread_ids[i] = i;
        pthread_create(&threads[i], NULL, pfunction, (void *)&thread_ids[i]);
    }

    for (int i = 0; i < T_PTHREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    // Destruir la barrera
    pthread_barrier_destroy(&barrier);

    MPI_Finalize();
    return 0;
}

/* Para tiempo de ejecucion */
double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

/* Pthreads */
void *pfunction(void *arg)
{
    int idW = *((int *)arg);
    int slice = N / T_PTHREADS;
    int ini = idW * slice;
    int lim = ini + slice;

    // Esperar en la barrera antes de imprimir
    pthread_barrier_wait(&barrier);

    printf("Thread %d (Proceso MPI %d): Porción del arreglo [%d, %d):\n", idW, idW_MPI, ini, lim);
    for (int i = ini; i < lim && i < N; i++)
    {
        printf("%d ", array[i]);
    }
    printf("\n");

    pthread_exit(NULL);
}