#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <string.h>

void calcularFuerzasEntreBloques(cuerpo_t *bloqueA, cuerpo_t *bloqueB, int sizeA, int sizeB,
                                 double *fuerzaX, double *fuerzaY, double *fuerzaZ, int offsetA)
{
    int i, j;
    double dif_X, dif_Y, dif_Z;
    double distancia, F;

    for (i = 0; i < sizeA; i++)
    {
        for (j = 0; j < sizeB; j++)
        {
            // Si están en la misma posición, ignorar
            if ((bloqueA[i].px == bloqueB[j].px) &&
                (bloqueA[i].py == bloqueB[j].py) &&
                (bloqueA[i].pz == bloqueB[j].pz))
                continue;

            dif_X = bloqueB[j].px - bloqueA[i].px;
            dif_Y = bloqueB[j].py - bloqueA[i].py;
            dif_Z = bloqueB[j].pz - bloqueA[i].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);
            F = (G * bloqueA[i].masa * bloqueB[j].masa) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            fuerzaX[offsetA + i] += dif_X;
            fuerzaY[offsetA + i] += dif_Y;
            fuerzaZ[offsetA + i] += dif_Z;

            fuerzaX[offsetA + sizeA + j] -= dif_X;
            fuerzaY[offsetA + sizeA + j] -= dif_Y;
            fuerzaZ[offsetA + sizeA + j] -= dif_Z;
        }
    }
}

void calcularFuerzasInternas(cuerpo_t *bloque, int size,
                             double *fuerzaX, double *fuerzaY, double *fuerzaZ,
                             int offset)
{
    int i, j;
    double dif_X, dif_Y, dif_Z;
    double distancia, F;

    for (i = 0; i < size - 1; i++)
    {
        for (j = i + 1; j < size; j++)
        {
            if ((bloque[i].px == bloque[j].px) &&
                (bloque[i].py == bloque[j].py) &&
                (bloque[i].pz == bloque[j].pz))
                continue;

            dif_X = bloque[j].px - bloque[i].px;
            dif_Y = bloque[j].py - bloque[i].py;
            dif_Z = bloque[j].pz - bloque[i].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);
            F = (G * bloque[i].masa * bloque[j].masa) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            fuerzaX[offset + i] += dif_X;
            fuerzaY[offset + i] += dif_Y;
            fuerzaZ[offset + i] += dif_Z;

            fuerzaX[offset + j] -= dif_X;
            fuerzaY[offset + j] -= dif_Y;
            fuerzaZ[offset + j] -= dif_Z;
        }
    }
}

void moverCuerpos(cuerpo_t *cuerpos, int N, int dt,
                  double *fuerzaX, double *fuerzaY, double *fuerzaZ,
                  int offset)
{
    for (int i = 0; i < N; i++)
    {
        fuerzaX[offset + i] *= 1.0 / cuerpos[i].masa;
        fuerzaY[offset + i] *= 1.0 / cuerpos[i].masa;
        // fuerzaZ[offset + i] *= 1.0 / cuerpos[i].masa;

        cuerpos[i].vx += fuerzaX[offset + i] * dt;
        cuerpos[i].vy += fuerzaY[offset + i] * dt;
        // cuerpos[i].vz += fuerzaZ[offset + i] * dt;

        cuerpos[i].px += cuerpos[i].vx * dt;
        cuerpos[i].py += cuerpos[i].vy * dt;
        // cuerpos[i].pz += cuerpos[i].vz * dt;

        // Reiniciar fuerzas (opcional si se hace fuera de esta función)
        fuerzaX[offset + i] = 0.0;
        fuerzaY[offset + i] = 0.0;
        fuerzaZ[offset + i] = 0.0;
    }
}


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

void calcularFuerzas(cuerpo_t *cuerpos, int N, int dt);

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
    int block_size = N / 2;

    // Inicialización y memoria
    cuerpo_t *all_cuerpos = malloc(sizeof(cuerpo_t) * N);
    cuerpo_t *block = all_cuerpos;
    fuerza_totalX = malloc(sizeof(double) * N);
    fuerza_totalY = malloc(sizeof(double) * N);
    fuerza_totalZ = malloc(sizeof(double) * N);
    inicializarCuerpos(all_cuerpos, N);

    // Enviar todos los cuerpos al Worker
    MPI_Bcast(all_cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++)
    {
        // Limpiar fuerzas
        memset(fuerza_totalX, 0, sizeof(double) * N);
        memset(fuerza_totalY, 0, sizeof(double) * N);
        memset(fuerza_totalZ, 0, sizeof(double) * N);

        // Calcular fuerzas entre bloque 0 (coordinador) y bloque 1 (worker)
        calcularFuerzasEntreBloques(block, all_cuerpos + block_size, block_size, block_size,
                                    fuerza_totalX, fuerza_totalY, fuerza_totalZ, 0);

        // Recibir fuerzas internas del Worker
        MPI_Recv(fuerza_totalX + block_size, block_size, MPI_DOUBLE, 1, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(fuerza_totalY + block_size, block_size, MPI_DOUBLE, 1, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(fuerza_totalZ + block_size, block_size, MPI_DOUBLE, 1, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Calcular fuerzas internas propias
        calcularFuerzasInternas(block, block_size, fuerza_totalX, fuerza_totalY, fuerza_totalZ, 0);

        // Mover cuerpos propios
        moverCuerpos(block, block_size, dt, fuerza_totalX, fuerza_totalY, fuerza_totalZ);

        // Recibir posiciones del Worker
        MPI_Recv(all_cuerpos + block_size, block_size * sizeof(cuerpo_t), MPI_BYTE, 1, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    tFin = dwalltime();
    tTotal = tFin - tIni;
    printf("Tiempo en segundos: %f\n", tTotal);

    for (int i = 0; i < N; i++)
    {
        printf("%f\n%f\n%f\n", all_cuerpos[i].px, all_cuerpos[i].py, all_cuerpos[i].pz);
    }

    finalizar();
}

void Worker(void)
{
    int block_size = N / 2;

    // Solo trabajamos con la mitad superior
    cuerpo_t *all_cuerpos = malloc(sizeof(cuerpo_t) * N);
    cuerpo_t *block = all_cuerpos + block_size;
    fuerza_totalX = malloc(sizeof(double) * N);
    fuerza_totalY = malloc(sizeof(double) * N);
    fuerza_totalZ = malloc(sizeof(double) * N);

    // Recibir todos los cuerpos
    MPI_Bcast(all_cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int paso = 0; paso < pasos; paso++)
    {
        memset(fuerza_totalX, 0, sizeof(double) * N);
        memset(fuerza_totalY, 0, sizeof(double) * N);
        memset(fuerza_totalZ, 0, sizeof(double) * N);

        // Calcular fuerzas internas de su propio bloque
        calcularFuerzasInternas(block, block_size,
                                fuerza_totalX + block_size,
                                fuerza_totalY + block_size,
                                fuerza_totalZ + block_size,
                                block_size);

        // Enviar fuerzas al Coordinador
        MPI_Send(fuerza_totalX + block_size, block_size, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
        MPI_Send(fuerza_totalY + block_size, block_size, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD);
        MPI_Send(fuerza_totalZ + block_size, block_size, MPI_DOUBLE, 0, 102, MPI_COMM_WORLD);

        // Mover cuerpos propios
        moverCuerpos(block, block_size, dt,
                     fuerza_totalX + block_size,
                     fuerza_totalY + block_size,
                     fuerza_totalZ + block_size);

        // Enviar cuerpos actualizados al Coordinador
        MPI_Send(block, block_size * sizeof(cuerpo_t), MPI_BYTE, 0, 200, MPI_COMM_WORLD);
    }
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