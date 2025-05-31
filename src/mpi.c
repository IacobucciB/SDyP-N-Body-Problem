#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <mpi.h>

double dwalltime();
double tIni, tFin, tTotal;

//
// Constantes para Algoritmo de gravitacion
//
#define PI (3.141592653589793)
#define M_PI (3.14159265358979323846)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 // Hidrogeno molecular

//
// Estructuras y variables para Algoritmo de gravitacion
//
typedef struct cuerpo cuerpo_t;
struct cuerpo
{
    float masa;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float r;
    float g;
    float b;
    int cuerpo;
};

float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
// Add vectors to store last positions
float *lastPositionX, *lastPositionY, *lastPositionZ;
float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

cuerpo_t *cuerpos;
int delta_tiempo = 1.0f; // Intervalo de tiempo, longitud de un paso
int pasos;
int N;

int rank;
int T;
int dt;

void calcularFuerzas(int ini, int lim);
void moverCuerpos(int ini, int lim);
void inicializarEstrella(cuerpo_t *cuerpo, int i, double n);
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n);
void inicializarH2(cuerpo_t *cuerpo, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos, int N);
void finalizar(void);

pthread_barrier_t barrera;
void *pfunction(void *arg);

void mpiworker();

int main(int argc, char const *argv[])
{
    if (argc < 5)
    {
        printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos> <threads>\n", argv[0]);
        return -1;
    }

    int id;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    int cantidadDeProcesos;
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Comm_size(MPI_COMM_WORLD,&cantidadDeProcesos);

    mpiworker();


    MPI_Finalize();
    return 0;
}

void mpiworker(int idW)
{
    int slice = N / T;
    int ini = idW * slice;
    int lim = ini + slice;

    // numero maximo de otros cuerpos en mensajes
    int tempSize = N - slice;
}