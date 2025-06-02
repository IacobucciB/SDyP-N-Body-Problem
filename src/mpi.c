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

/* Argumentos */
int N;        // Número de cuerpos
int dt;       // Intervalo de tiempo, longitud de un paso
int pasos;    // Número de pasos a simular
int T_MPI;   // Número total de procesos MPI
int T_PTHREADS; // Número de threads por proceso

/* MAIN */
int main(int argc, char *argv[])
{
    // MPI
    int idW_MPI;
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
    int T_PTHREADS = atoi(argv[4]);

    if (N <= 0 || dt <= 0 || pasos <= 0 || T_PTHREADS <= 0) {
        fprintf(stderr, "Error: All arguments must be positive integers.\n");
        MPI_Finalize();
        return -1;
    }

    printf("N: %d, DT: %d, PASOS: %d, T_PTHREADS: %d\n", N, dt, pasos, T_PTHREADS);

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