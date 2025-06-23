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

void calcularFuerzas(cuerpo_t *all_cuerpos, double *local_fuerza_X, double *local_fuerza_Y, double *local_fuerza_Z, int N_total, int current_process_ini, int current_process_lim);

void moverCuerpos(cuerpo_t *cuerpos, int N, int dt);

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
    // These will hold the *final, summed* forces
    fuerza_totalX = malloc(sizeof(double) * N);
    fuerza_totalY = malloc(sizeof(double) * N);
    fuerza_totalZ = malloc(sizeof(double) * N);

    // Each process needs a temporary array to store the forces it calculates
    // for its assigned bodies.
    double *temp_fuerzaX = malloc(sizeof(double) * N);
    double *temp_fuerzaY = malloc(sizeof(double) * N);
    double *temp_fuerzaZ = malloc(sizeof(double) * N);

    inicializarCuerpos(cuerpos, N);
    MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes have 'cuerpos' before timing

    tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++)
    {
        // Before each step, clear the *final* force arrays on the root
        memset(fuerza_totalX, 0, sizeof(double) * N);
        memset(fuerza_totalY, 0, sizeof(double) * N);
        memset(fuerza_totalZ, 0, sizeof(double) * N);

        // All processes (including Coordinator) calculate forces for their assigned range
        calcularFuerzas(cuerpos, temp_fuerzaX, temp_fuerzaY, temp_fuerzaZ, N, ini_MPI, lim_MPI);

        // Sum up the force contributions from all processes to the root's fuerza_total arrays
        // The root (0) receives the sum of all temp_fuerza arrays into its fuerza_total arrays.
        MPI_Reduce(temp_fuerzaX, fuerza_totalX, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(temp_fuerzaY, fuerza_totalY, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(temp_fuerzaZ, fuerza_totalZ, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        // Only the Coordinator moves the bodies as it has the final, summed forces
        moverCuerpos(cuerpos, N, dt);

        // Broadcast the updated body positions to all workers for the next step
        MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize for next step
    }

    tFin = dwalltime();
    tTotal = tFin - tIni;
    printf("Tiempo en segundos: %f\n", tTotal);
    for (int i = 0; i < N; i++)
    {
        printf("%f\n%f\n%f\n", cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
    }

    free(temp_fuerzaX);
    free(temp_fuerzaY);
    free(temp_fuerzaZ);
    finalizar(); // Frees cuerpos, fuerza_totalX/Y/Z
}

void Worker(void)
{
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    // These will hold the local force contributions before reduction
    fuerza_totalX = malloc(sizeof(double) * N); // Still need for moverCuerpos logic later if not refactored
    fuerza_totalY = malloc(sizeof(double) * N);
    fuerza_totalZ = malloc(sizeof(double) * N);

    double *temp_fuerzaX = malloc(sizeof(double) * N);
    double *temp_fuerzaY = malloc(sizeof(double) * N);
    double *temp_fuerzaZ = malloc(sizeof(double) * N);

    MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int paso = 0; paso < pasos; paso++)
    {
        // Each worker calculates forces for its assigned range
        calcularFuerzas(cuerpos, temp_fuerzaX, temp_fuerzaY, temp_fuerzaZ, N, ini_MPI, lim_MPI);

        // All workers send their contributions to the root (rank 0) for summation
        MPI_Reduce(temp_fuerzaX, NULL, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // NULL for non-root receive buffer
        MPI_Reduce(temp_fuerzaY, NULL, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(temp_fuerzaZ, NULL, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        // Workers receive the updated body positions from the Coordinator
        MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize for next step
    }

    free(temp_fuerzaX);
    free(temp_fuerzaY);
    free(temp_fuerzaZ);
    finalizar(); // Frees cuerpos, fuerza_totalX/Y/Z (allocated in this process)
}

void calcularFuerzas(cuerpo_t *all_cuerpos, double *local_fuerza_X, double *local_fuerza_Y, double *local_fuerza_Z, int N_total, int current_process_ini, int current_process_lim)
{
    double dif_X, dif_Y, dif_Z;
    double distancia_sq, distancia;
    double F;

    // Initialize local force accumulators to zero for the assigned range
    for (int i = current_process_ini; i < current_process_lim; i++) {
        local_fuerza_X[i] = 0.0;
        local_fuerza_Y[i] = 0.0;
        local_fuerza_Z[i] = 0.0;
    }

    for (int i = current_process_ini; i < current_process_lim; i++)
    {
        for (int j = 0; j < N_total; j++)
        {
            if (i == j)
                continue;

            dif_X = all_cuerpos[j].px - all_cuerpos[i].px;
            dif_Y = all_cuerpos[j].py - all_cuerpos[i].py;
            dif_Z = all_cuerpos[j].pz - all_cuerpos[i].pz;

            distancia_sq = dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z;

            if (distancia_sq == 0.0)
                continue;

            distancia = sqrt(distancia_sq);

            F = (G * all_cuerpos[i].masa * all_cuerpos[j].masa) / (distancia_sq * distancia);

            local_fuerza_X[i] += F * (dif_X / distancia);
            local_fuerza_Y[i] += F * (dif_Y / distancia);
            local_fuerza_Z[i] += F * (dif_Z / distancia);
        }
    }
}

void moverCuerpos(cuerpo_t *cuerpos, int N, int dt)
{
	int cuerpo;
	for (cuerpo = 0; cuerpo < N; cuerpo++)
	{

		fuerza_totalX[cuerpo] *= 1 / cuerpos[cuerpo].masa;
		fuerza_totalY[cuerpo] *= 1 / cuerpos[cuerpo].masa;
		// fuerza_totalZ[cuerpo] *= 1/cuerpos[cuerpo].masa;

		cuerpos[cuerpo].vx += fuerza_totalX[cuerpo] * dt;
		cuerpos[cuerpo].vy += fuerza_totalY[cuerpo] * dt;
		// cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*dt;

		cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;
		cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;
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