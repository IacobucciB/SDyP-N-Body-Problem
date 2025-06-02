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

void calcularFuerzas(int ini, int lim, float *local_fuerzaX, float *local_fuerzaY, float *local_fuerzaZ)
{
    int cuerpo1, cuerpo2;
    float dif_X, dif_Y, dif_Z;
    float distancia;
    float F;

    for (cuerpo1 = ini; cuerpo1 < lim; cuerpo1++)
    {
        for (cuerpo2 = 0; cuerpo2 < N; cuerpo2++)
        {
            if (cuerpo1 == cuerpo2) continue;
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

            // Store forces in local arrays
            local_fuerzaX[cuerpo1] += dif_X;
            local_fuerzaY[cuerpo1] += dif_Y;
            local_fuerzaZ[cuerpo1] += dif_Z;
        }
    }
}

void moverCuerpos(int ini, int lim)
{
    int cuerpo;
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

    srand(time(NULL));

    for (cuerpo = 0; cuerpo < N; cuerpo++)
    {

        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;

        cuerpos[cuerpo].cuerpo = (rand() % 3);

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
    if (cuerpos) free(cuerpos);
    if (fuerza_totalX) free(fuerza_totalX);
    if (fuerza_totalY) free(fuerza_totalY);
    if (fuerza_totalZ) free(fuerza_totalZ);
    if (lastPositionX) free(lastPositionX);
    if (lastPositionY) free(lastPositionY);
    if (lastPositionZ) free(lastPositionZ);
}

pthread_barrier_t barrera;

void *pfunction(void *arg)
{
    int idW = *((int *)arg);
    int paso;
    int slice = N / T;
    int ini = idW * slice;
    int lim = ini + slice;

    // Allocate full-size local force arrays using malloc
    float *local_fuerzaX = (float *)malloc(N * sizeof(float));
    float *local_fuerzaY = (float *)malloc(N * sizeof(float));
    float *local_fuerzaZ = (float *)malloc(N * sizeof(float));

    if (!local_fuerzaX || !local_fuerzaY || !local_fuerzaZ) {
        fprintf(stderr, "Error al asignar memoria para las matrices locales de fuerzas.\n");
        if (local_fuerzaX) free(local_fuerzaX);
        if (local_fuerzaY) free(local_fuerzaY);
        if (local_fuerzaZ) free(local_fuerzaZ);
        pthread_exit(NULL);
    }

    // Initialize arrays to zero
    for (int i = 0; i < N; i++) {
        local_fuerzaX[i] = 0.0f;
        local_fuerzaY[i] = 0.0f;
        local_fuerzaZ[i] = 0.0f;
    }

    for (paso = 0; paso < pasos; paso++)
    {
        calcularFuerzas(ini, lim, local_fuerzaX, local_fuerzaY, local_fuerzaZ);
        pthread_barrier_wait(&barrera);

        // Combine forces from all threads
        for (int i = ini; i < lim; i++)
        {
            fuerza_totalX[i] = local_fuerzaX[i];
            fuerza_totalY[i] = local_fuerzaY[i];
            fuerza_totalZ[i] = local_fuerzaZ[i];
        }

        pthread_barrier_wait(&barrera);
        moverCuerpos(ini, lim);

        // Reset local forces arrays using for loop
        for (int i = 0; i < N; i++) {
            local_fuerzaX[i] = 0.0f;
            local_fuerzaY[i] = 0.0f;
            local_fuerzaZ[i] = 0.0f;
        }
    }

    free(local_fuerzaX);
    free(local_fuerzaY);
    free(local_fuerzaZ);

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
    fuerza_totalX = (float *)malloc(sizeof(float) * N);
    fuerza_totalY = (float *)malloc(sizeof(float) * N);
    fuerza_totalZ = (float *)malloc(sizeof(float) * N);
    // Allocate memory for last position vectors
    lastPositionX = (float *)malloc(sizeof(float) * N);
    lastPositionY = (float *)malloc(sizeof(float) * N);
    lastPositionZ = (float *)malloc(sizeof(float) * N);


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
    // Save last positions after simulation
    for (int i = 0; i < N; i++)
    {
        lastPositionX[i] = cuerpos[i].px;
        lastPositionY[i] = cuerpos[i].py;
        lastPositionZ[i] = cuerpos[i].pz;
    }

    tTotal = tFin - tIni;

    //Print last positions of all bodies
    printf("\n=== Last Positions of Bodies ===\n");
    printf("%-6s %-15s %-15s %-15s\n", "ID", "X", "Y", "Z");
    for (int i = 0; i < N; i++)
    {
        printf("%-6d %-15.6f %-15.6f %-15.6f\n", i, lastPositionX[i], lastPositionY[i], lastPositionZ[i]);
    }

    printf("Tiempo en segundos: %f\n", tTotal);

    pthread_barrier_destroy(&barrera);
    finalizar();
    return 0;
}
