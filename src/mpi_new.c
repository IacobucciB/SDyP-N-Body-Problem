#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <string.h>

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

int T_PTHREADS;

void *pcalcularFuerzas(void *arg);
void *pcalcularFuerzasB0(void *arg);
void *pcalcularFuerzasB1(void *arg);
void *pcalcularFuerzasCruzadas(void *arg);
void *pmoverCuerpos(void *arg);

void calcularFuerzas(cuerpo_t *cuerpos, int N, int dt);

// void moverCuerpos(cuerpo_t *cuerpos, int N, int dt);
void moverCuerpos(cuerpo_t *cuerpos, int inicio, int fin, int dt);

int idW_MPI;
int T_MPI;
int blockSize;
int ini_MPI;
int lim_MPI;
int tempSize;
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

    blockSize = N / T_MPI;
    ini_MPI = idW_MPI * blockSize;
    lim_MPI = ini_MPI + blockSize;
    tempSize = N - blockSize;

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
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (double *)malloc(sizeof(double) * N * T_PTHREADS);
    fuerza_totalY = (double *)malloc(sizeof(double) * N * T_PTHREADS);
    fuerza_totalZ = (double *)malloc(sizeof(double) * N * T_PTHREADS);

    inicializarCuerpos(cuerpos, N);

    MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // int T_PTHREADS;
    pthread_t threads[T_PTHREADS];
    int thread_ids[T_PTHREADS];
    int mid = N / 2;
    int resto = N - mid;
    tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++)
    {
        for (int i = 0; i < N * T_PTHREADS; i++)
        {
            fuerza_totalX[i] = 0.0;
            fuerza_totalY[i] = 0.0;
            fuerza_totalZ[i] = 0.0;
        }
        // calcularFuerzas(cuerpos, N, dt);
        for (int i = 0; i < T_PTHREADS; i++)
        {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, pcalcularFuerzasB0, &thread_ids[i]);
        }
        for (int i = 0; i < T_PTHREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        // MPI_Bcast(fuerza_totalX, N * T_PTHREADS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // MPI_Bcast(fuerza_totalY, N * T_PTHREADS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // MPI_Bcast(fuerza_totalZ, N * T_PTHREADS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // moverCuerpos(cuerpos, 0, mid, dt);
        MPI_Recv(&cuerpos[mid], resto * sizeof(cuerpo_t), MPI_BYTE, 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < T_PTHREADS; i++)
        {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, pcalcularFuerzasCruzadas, &thread_ids[i]);
        }
        for (int i = 0; i < T_PTHREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        MPI_Send(fuerza_totalX[mid], resto * sizeof(double), MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
        MPI_Send(fuerza_totalY[mid], resto * sizeof(double), MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
        MPI_Send(fuerza_totalZ[mid], resto * sizeof(double), MPI_DOUBLE, 1, 3, MPI_COMM_WORLD);

        for (int i = 0; i < T_PTHREADS; i++)
        {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, pmoverCuerpos, &thread_ids[i]);
        }
        for (int i = 0; i < T_PTHREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
       // MPI_Recv(&cuerpos[mid], resto * sizeof(cuerpo_t), MPI_BYTE, 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    tFin = dwalltime();
    tTotal = tFin - tIni;
    printf("%f\n", tTotal);

    for (int i = 0; i < N; i++)
    {
        printf("%f\n%f\n%f\n", cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
    }

    finalizar();
}

void Worker(void)
{
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (double *)malloc(sizeof(double) * N * T_PTHREADS);
    fuerza_totalY = (double *)malloc(sizeof(double) * N * T_PTHREADS);
    fuerza_totalZ = (double *)malloc(sizeof(double) * N * T_PTHREADS);
    MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    int mid = N / 2;
    int resto = N - mid;
    pthread_t threads[T_PTHREADS];
    int thread_ids[T_PTHREADS];
    for (int paso = 0; paso < pasos; paso++)
    {
        for (int i = 0; i < N * T_PTHREADS; i++)
        {
            fuerza_totalX[i] = 0.0;
            fuerza_totalY[i] = 0.0;
            fuerza_totalZ[i] = 0.0;
        }
        MPI_Send(&cuerpos[mid], resto * sizeof(cuerpo_t), MPI_BYTE, 0, 4, MPI_COMM_WORLD);

        for (int i = 0; i < T_PTHREADS; i++)
        {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, pcalcularFuerzasB1, &thread_ids[i]);
        }
        for (int i = 0; i < T_PTHREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        double *tfX = (double *)malloc(resto * sizeof(double));
        double *tfY = (double *)malloc(resto * sizeof(double));
        double *tfZ = (double *)malloc(resto * sizeof(double));
        MPI_Recv(tfX, resto * sizeof(double), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(tfY, resto * sizeof(double), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(tfZ, resto * sizeof(double), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < resto; i++)
        {
            fuerza_totalX[idW_MPI * N + mid + i] += tfX[i];
            fuerza_totalY[idW_MPI * N + mid + i] += tfY[i];
            fuerza_totalZ[idW_MPI * N + mid + i] += tfZ[i];
        }

        // MPI_Bcast(fuerza_totalX, N * T_PTHREADS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // MPI_Bcast(fuerza_totalY, N * T_PTHREADS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // MPI_Bcast(fuerza_totalZ, N * T_PTHREADS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // moverCuerpos(cuerpos, mid, N, dt);
        for (int i = 0; i < T_PTHREADS; i++)
        {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, pmoverCuerpos, &thread_ids[i]);
        }
        for (int i = 0; i < T_PTHREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        // MPI_Send(&cuerpos[mid], resto * sizeof(cuerpo_t), MPI_BYTE, 0, 4, MPI_COMM_WORLD);
    }
    finalizar();
}

void calcularFuerzas(cuerpo_t *cuerpos, int N, int dt)
{
    int cuerpo1, cuerpo2;
    double dif_X, dif_Y, dif_Z;
    double distancia;
    double F;

    for (cuerpo1 = 0; cuerpo1 < N - 1; cuerpo1++)
    {
        for (cuerpo2 = cuerpo1 + 1; cuerpo2 < N; cuerpo2++)
        {
            if ((cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) && (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) && (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

            dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
            dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
            dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            F = (G * cuerpos[cuerpo1].masa * cuerpos[cuerpo2].masa) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            fuerza_totalX[cuerpo1] += dif_X;
            fuerza_totalY[cuerpo1] += dif_Y;
            fuerza_totalZ[cuerpo1] += dif_Z;

            fuerza_totalX[cuerpo2] -= dif_X;
            fuerza_totalY[cuerpo2] -= dif_Y;
            fuerza_totalZ[cuerpo2] -= dif_Z;
        }
    }
}

void *pcalcularFuerzas(void *arg)
{
    int idW = *((int *)arg);
    int slice = N / T_PTHREADS;
    int ini = idW * slice;
    int lim = ini + slice;

    int cuerpo1, cuerpo2;
    double dif_X, dif_Y, dif_Z;
    double distancia;
    double F;

    for (cuerpo1 = ini; cuerpo1 < lim; cuerpo1++)
    {
        for (cuerpo2 = cuerpo1 + 1; cuerpo2 < N; cuerpo2++)
        {
            if ((cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) &&
                (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) &&
                (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

            dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
            dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
            dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);
            F = (G * cuerpos[cuerpo1].masa * cuerpos[cuerpo2].masa) / (distancia * distancia);

            fuerza_totalX[idW * N + cuerpo1] += dif_X * F;
            fuerza_totalY[idW * N + cuerpo1] += dif_Y * F;
            fuerza_totalZ[idW * N + cuerpo1] += dif_Z * F;

            fuerza_totalX[idW * N + cuerpo2] -= dif_X * F;
            fuerza_totalY[idW * N + cuerpo2] -= dif_Y * F;
            fuerza_totalZ[idW * N + cuerpo2] -= dif_Z * F;
        }
    }

    pthread_exit(NULL);
}

void *pcalcularFuerzasCruzadas(void *arg)
{
    int tid = *((int *)arg);
    int slice = (N / 2) / T_PTHREADS;  // Solo iteramos B0
    int ini = tid * slice;
    int fin = ini + slice;

    int i, j;
    double dx, dy, dz, dist, F;

    for (i = ini; i < fin; i++)
    {
        for (j = N / 2; j < N; j++) // B1
        {
            dx = cuerpos[j].px - cuerpos[i].px;
            dy = cuerpos[j].py - cuerpos[i].py;
            dz = cuerpos[j].pz - cuerpos[i].pz;

            dist = sqrt(dx*dx + dy*dy + dz*dz);
            F = (G * cuerpos[i].masa * cuerpos[j].masa) / (dist * dist);

            dx *= F;
            dy *= F;
            dz *= F;

            fuerza_totalX[tid * N + i] += dx;
            fuerza_totalY[tid * N + i] += dy;
            fuerza_totalZ[tid * N + i] += dz;

            fuerza_totalX[tid * N + j] -= dx;
            fuerza_totalY[tid * N + j] -= dy;
            fuerza_totalZ[tid * N + j] -= dz;
        }
    }

    pthread_exit(NULL);
}

void *pcalcularFuerzasB0(void *arg)
{
    int tid = *((int *)arg);
    int mid = N / 2;
    int slice = mid / T_PTHREADS;
    int ini = tid * slice;
    int fin = ini + slice;

    int i, j;
    double dx, dy, dz, dist, F;

    for (i = ini; i < fin; i++) {
        for (j = i + 1; j < mid; j++) {  // Solo hasta mid-1
            dx = cuerpos[j].px - cuerpos[i].px;
            dy = cuerpos[j].py - cuerpos[i].py;
            dz = cuerpos[j].pz - cuerpos[i].pz;

            dist = sqrt(dx * dx + dy * dy + dz * dz);
            F = (G * cuerpos[i].masa * cuerpos[j].masa) / (dist * dist);

            dx *= F;
            dy *= F;
            dz *= F;

            fuerza_totalX[tid * N + i] += dx;
            fuerza_totalY[tid * N + i] += dy;
            fuerza_totalZ[tid * N + i] += dz;

            fuerza_totalX[tid * N + j] -= dx;
            fuerza_totalY[tid * N + j] -= dy;
            fuerza_totalZ[tid * N + j] -= dz;
        }
    }

    pthread_exit(NULL);
}

void *pcalcularFuerzasB1(void *arg)
{
    int tid = *((int *)arg);
    int mid = N / 2;
    int resto = N - mid;
    int slice = resto / T_PTHREADS;
    int ini = mid + tid * slice;
    int fin = ini + slice;

    int i, j;
    double dx, dy, dz, dist, F;

    for (i = ini; i < fin; i++) {
        for (j = i + 1; j < N; j++) {
            dx = cuerpos[j].px - cuerpos[i].px;
            dy = cuerpos[j].py - cuerpos[i].py;
            dz = cuerpos[j].pz - cuerpos[i].pz;

            dist = sqrt(dx * dx + dy * dy + dz * dz);
            F = (G * cuerpos[i].masa * cuerpos[j].masa) / (dist * dist);

            dx *= F;
            dy *= F;
            dz *= F;

            fuerza_totalX[tid * N + i] += dx;
            fuerza_totalY[tid * N + i] += dy;
            fuerza_totalZ[tid * N + i] += dz;

            fuerza_totalX[tid * N + j] -= dx;
            fuerza_totalY[tid * N + j] -= dy;
            fuerza_totalZ[tid * N + j] -= dz;
        }
    }

    pthread_exit(NULL);
}


// void moverCuerpos(cuerpo_t *cuerpos, int N, int dt)
void moverCuerpos(cuerpo_t *cuerpos, int inicio, int fin, int dt)
{
    int cuerpo;
    for (cuerpo = inicio; cuerpo < fin; cuerpo++)
    {
        double fx = 0.0, fy = 0.0, fz = 0.0;

        for (int t = 0; t < T_PTHREADS; t++)
        {
            fx += fuerza_totalX[t * N + cuerpo];
            fy += fuerza_totalY[t * N + cuerpo];
            fz += fuerza_totalZ[t * N + cuerpo];
        }

        fx /= cuerpos[cuerpo].masa;
        fy /= cuerpos[cuerpo].masa;
        // fz /= cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fx * dt;
        cuerpos[cuerpo].vy += fy * dt;
        // cuerpos[cuerpo].vz += fz * delta_tiempo;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;
        // cuerpos[cuerpo].pz += cuerpos[cuerpo].vz * delta_tiempo;

        for (int t = 0; t < T_PTHREADS; t++)
        {
            fuerza_totalX[t * N + cuerpo] = 0.0;
            fuerza_totalY[t * N + cuerpo] = 0.0;
            fuerza_totalZ[t * N + cuerpo] = 0.0;
        }
    }
}

void *pmoverCuerpos(void *arg)
{
    int idW = *((int *)arg);
    int slice = N / T_PTHREADS;
    int ini = idW * slice;
    int lim = ini + slice;

    if (idW_MPI == 0 && ini == 0)
    {
        ini = 0;
        lim = N / 2;
    }
    else if (idW_MPI == 1 && ini > 0)
    {
        ini = N / 2;
        lim = N;
    }

    for (int cuerpo = ini; cuerpo < lim; cuerpo++)
    {
        double fx = 0.0, fy = 0.0, fz = 0.0;

        for (int t = 0; t < T_PTHREADS; t++)
        {
            fx += fuerza_totalX[t * N + cuerpo];
            fy += fuerza_totalY[t * N + cuerpo];
            fz += fuerza_totalZ[t * N + cuerpo];
        }

        fx /= cuerpos[cuerpo].masa;
        fy /= cuerpos[cuerpo].masa;
        // fz /= cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fx * dt;
        cuerpos[cuerpo].vy += fy * dt;
        // cuerpos[cuerpo].vz += fz * delta_tiempo;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;
        // cuerpos[cuerpo].pz += cuerpos[cuerpo].vz * delta_tiempo;

        for (int t = 0; t < T_PTHREADS; t++)
        {
            fuerza_totalX[t * N + cuerpo] = 0.0;
            fuerza_totalY[t * N + cuerpo] = 0.0;
            fuerza_totalZ[t * N + cuerpo] = 0.0;
        }
    }

    pthread_exit(NULL);
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

double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

void finalizar(void)
{
    free(cuerpos);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
}