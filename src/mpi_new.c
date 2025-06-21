#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

double dwalltime();
double tIni, tFin, tTotal;

#define PI (3.141592653589793)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2

typedef struct cuerpo cuerpo_t;
struct cuerpo
{
    double masa;
    double px;
    double py;
    double pz;
    double vx;
    double vy;
    double vz;
    double r;
    double g;
    double b;
    int cuerpo;
};

double *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
double toroide_alfa;
double toroide_theta;
double toroide_incremento;
double toroide_lado;
double toroide_r;
double toroide_R;

cuerpo_t *cuerpos;
double dt = 1.0f;
int pasos;
int N;

void inicializarEstrella(cuerpo_t *cuerpo, int i, double n);
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n);
void inicializarH2(cuerpo_t *cuerpo, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos, int N);

void finalizar(void);

int idW_MPI;
int T_MPI;
int T_PTHREADS;

void *calcularFuerzas(void *arg);
void *moverCuerpos(void *arg);

void Coordinator(void);
void Worker(void);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &idW_MPI);
    MPI_Comm_size(MPI_COMM_WORLD, &T_MPI);

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

    if (idW_MPI == 0)
    {
        Coordinator();
    }
    else
    {
        Worker();
    }

    MPI_Finalize();
    return 0;
}

void Coordinator(void)
{
    int slice_MPI = N / T_MPI;
    int ini_MPI = idW_MPI * slice_MPI;
    int lim_MPI = ini_MPI + slice_MPI;

    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (double *)malloc(sizeof(double) * N);
    fuerza_totalY = (double *)malloc(sizeof(double) * N);
    fuerza_totalZ = (double *)malloc(sizeof(double) * N);
    inicializarCuerpos(cuerpos, N);

    pthread_t threads[T_PTHREADS];
    int thread_ids[T_PTHREADS];

    MPI_Barrier(MPI_COMM_WORLD);

    tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++)
    {
    }

    tFin = dwalltime();
    tTotal = tFin - tIni;
    printf("Tiempo en segundos: %f\n", tTotal);
    for (int i = 0; i < N; i++)
    {
        printf("%f\n%f\n%f\n", cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
    }

    finalizar();
}

void Worker(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    for (int paso = 0; paso < pasos; paso++)
    {
        
    }
    
}