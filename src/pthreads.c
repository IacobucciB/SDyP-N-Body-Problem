#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

//
// Para tiempo de ejecucion
//

double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}
double tIni, tFin, tTotal;

//
// Constantes para Algoritmo de gravitacion
//
#define PI (3.141592653589793)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 // Hidrogeno molecular

// ===============
// ===== CPU =====
// ===============

//
// Estructuras y variables para Algoritmo de gravitacion
//
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
int delta_tiempo = 1.0f; // Intervalo de tiempo, longitud de un paso
int pasos;
int N;
int T;

// void calcularFuerzas(int ini, int lim, double *local_fuerzaX, double *local_fuerzaY, double *local_fuerzaZ)
void calcularFuerzas(int idW)
{
    int cuerpo1, cuerpo2;
    double dif_X, dif_Y, dif_Z;
    double distancia;
    double F;

    int slice = N / T;
    int ini = idW * slice;
    int lim = ini + slice;

    for (cuerpo1 = ini; cuerpo1 < lim; cuerpo1++)
    {
        for (cuerpo2 = cuerpo1 + 1; cuerpo2 < N; cuerpo2++)
        {
            if (cuerpo1 == cuerpo2)
                continue;
            if ((cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) &&
                (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) &&
                (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

            dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
            dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
            dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);
            F = (G * cuerpos[cuerpo1].masa * cuerpos[cuerpo2].masa) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            /*
            fuerza_totalX[idW * N + cuerpo1] += dif_X;
            fuerza_totalY[idW * N + cuerpo1] += dif_Y;
            fuerza_totalZ[idW * N + cuerpo1] += dif_Z;
            fuerza_totalX[idW * N + cuerpo2] -= dif_X;
            fuerza_totalY[idW * N + cuerpo2] -= dif_Y;
            fuerza_totalZ[idW * N + cuerpo2] -= dif_Z;
            */
            
            fuerza_totalX[cuerpo1] += dif_X;
            fuerza_totalY[cuerpo1] += dif_Y;
            fuerza_totalZ[cuerpo1] += dif_Z;

            fuerza_totalX[cuerpo2] -= dif_X;
            fuerza_totalY[cuerpo2] -= dif_Y;
            fuerza_totalZ[cuerpo2] -= dif_Z;
        }
    }
}

void moverCuerpos(int idW)
{
    int cuerpo;
    int slice = N / T;
    int ini = idW * slice;
    int lim = ini + slice;

    for (cuerpo = ini; cuerpo < lim; cuerpo++)
    {

        fuerza_totalX[cuerpo] *= 1 / cuerpos[cuerpo].masa;
        fuerza_totalY[cuerpo] *= 1 / cuerpos[cuerpo].masa;
        // fuerza_totalZ[cuerpo] *= 1/cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo] * delta_tiempo;
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo] * delta_tiempo;
        // cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * delta_tiempo;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * delta_tiempo;
        // cuerpos[cuerpo].pz += cuerpos[cuerpo].vz *dt;

        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;
    }
}

void inicializarEstrella(cuerpo_t *cuerpo, int i, double n)
{

    cuerpo->masa = 0.001 * 8;

    if ((toroide_alfa + toroide_incremento) >= 2 * M_PI)
    {
        toroide_alfa = 0;
        toroide_theta += toroide_incremento;
    }
    else
    {
        toroide_alfa += toroide_incremento;
    }

    cuerpo->px = (toroide_R + toroide_r * cos(toroide_alfa)) * cos(toroide_theta);
    cuerpo->py = (toroide_R + toroide_r * cos(toroide_alfa)) * sin(toroide_theta);
    cuerpo->pz = toroide_r * sin(toroide_alfa);

    cuerpo->vx = 0.0;
    cuerpo->vy = 0.0;
    cuerpo->vz = 0.0;

    cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
    cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
    cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarPolvo(cuerpo_t *cuerpo, int i, double n)
{

    cuerpo->masa = 0.001 * 4;

    if ((toroide_alfa + toroide_incremento) >= 2 * M_PI)
    {
        toroide_alfa = 0;
        toroide_theta += toroide_incremento;
    }
    else
    {
        toroide_alfa += toroide_incremento;
    }

    cuerpo->px = (toroide_R + toroide_r * cos(toroide_alfa)) * cos(toroide_theta);
    cuerpo->py = (toroide_R + toroide_r * cos(toroide_alfa)) * sin(toroide_theta);
    cuerpo->pz = toroide_r * sin(toroide_alfa);

    cuerpo->vx = 0.0;
    cuerpo->vy = 0.0;
    cuerpo->vz = 0.0;

    cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
    cuerpo->g = 0.0; //(double )rand()/(RAND_MAX+1.0);
    cuerpo->b = 0.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarH2(cuerpo_t *cuerpo, int i, double n)
{

    cuerpo->masa = 0.001;

    if ((toroide_alfa + toroide_incremento) >= 2 * M_PI)
    {
        toroide_alfa = 0;
        toroide_theta += toroide_incremento;
    }
    else
    {
        toroide_alfa += toroide_incremento;
    }

    cuerpo->px = (toroide_R + toroide_r * cos(toroide_alfa)) * cos(toroide_theta);
    cuerpo->py = (toroide_R + toroide_r * cos(toroide_alfa)) * sin(toroide_theta);
    cuerpo->pz = toroide_r * sin(toroide_alfa);

    cuerpo->vx = 0.0;
    cuerpo->vy = 0.0;
    cuerpo->vz = 0.0;

    cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
    cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
    cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarCuerpos(cuerpo_t *cuerpos, int N)
{
    int cuerpo;
    double n = N;

    toroide_alfa = 0.0;
    toroide_theta = 0.0;
    toroide_lado = sqrt(N);
    toroide_incremento = 2 * M_PI / toroide_lado;
    toroide_r = 1.0;
    toroide_R = 2 * toroide_r;

    srand(0);

    for (cuerpo = 0; cuerpo < N; cuerpo++)
    {

        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;

        // cuerpos[cuerpo].cuerpo = (rand() % 3);
        cuerpos[cuerpo].cuerpo = cuerpo % 3;

        if (cuerpos[cuerpo].cuerpo == ESTRELLA)
        {
            inicializarEstrella(&cuerpos[cuerpo], cuerpo, n);
        }
        else if (cuerpos[cuerpo].cuerpo == POLVO)
        {
            inicializarPolvo(&cuerpos[cuerpo], cuerpo, n);
        }
        else if (cuerpos[cuerpo].cuerpo == H2)
        {
            inicializarH2(&cuerpos[cuerpo], cuerpo, n);
        }
    }

    cuerpos[0].masa = 2.0e2;
    cuerpos[0].px = 0.0;
    cuerpos[0].py = 0.0;
    cuerpos[0].pz = 0.0;
    cuerpos[0].vx = -0.000001;
    cuerpos[0].vy = -0.000001;
    cuerpos[0].vz = 0.0;

    cuerpos[1].masa = 1.0e1;
    cuerpos[1].px = -1.0;
    cuerpos[1].py = 0.0;
    cuerpos[1].pz = 0.0;
    cuerpos[1].vx = 0.0;
    cuerpos[1].vy = 0.0001;
    cuerpos[1].vz = 0.0;
}

void finalizar(void)
{
    free(cuerpos);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
}

pthread_barrier_t barrera;

void *pfunction(void *arg)
{
    int idW = *((int *)arg);
    int paso;
    int slice = N / T;
    int ini = idW * slice;
    int lim = ini + slice;

    for (paso = 0; paso < pasos; paso++)
    {
        calcularFuerzas(idW);
        pthread_barrier_wait(&barrera);
        moverCuerpos(idW);
        pthread_barrier_wait(&barrera);
    }
    pthread_exit(NULL);
}

int main(int argc, char const *argv[])
{
    if (argc < 5)
    {
        printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos> <threads>\n", argv[0]);
        return -1;
    }

    N = atoi(argv[1]);
    delta_tiempo = atof(argv[2]);
    pasos = atoi(argv[3]);
    T = atoi(argv[4]);

    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (double *)malloc(sizeof(double) * N);
    fuerza_totalY = (double *)malloc(sizeof(double) * N);
    fuerza_totalZ = (double *)malloc(sizeof(double) * N);

    inicializarCuerpos(cuerpos, N);

    pthread_t threads[T];
    int thread_ids[T];
    pthread_barrier_init(&barrera, NULL, T);

    tIni = dwalltime();

    for (int i = 0; i < T; i++)
    {
        thread_ids[i] = i;
        pthread_create(&threads[i], NULL, pfunction, &thread_ids[i]);
    }

    for (int i = 0; i < T; i++)
    {
        pthread_join(threads[i], NULL);
    }
    tFin = dwalltime();
    tTotal = tFin - tIni;

    printf("Tiempo en segundos: %f\n", tTotal);
    for (int i = 0; i < N; i++)
    {
        printf("%f\n%f\n%f\n", cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
    }

    pthread_barrier_destroy(&barrera);
    finalizar();
    return 0;
}
