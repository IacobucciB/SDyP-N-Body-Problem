#ifndef NBODY_H
#define NBODY_H

#include <pthread.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// #define PI (3.141592653589793)
#define M_PI (3.14159265358979323846264338327950288)
#define G 6.673e-3
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
double dt;
int pasos;
int N;
int T;

void calcularFuerzas(int idW);
void moverCuerpos(int idW);

void inicializarEstrella(cuerpo_t *cuerpo, int i, double n);
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n);
void inicializarH2(cuerpo_t *cuerpo, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos, int N);

void finalizar(void);

// pthread_t *threads;
// int *thread_ids;
pthread_barrier_t barrera;
void *pfunction(void *arg);

#endif // NBODY_H