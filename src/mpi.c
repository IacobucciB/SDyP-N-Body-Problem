#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <mpi.h>

double dwalltime();
double tIni, tFin, tTotal;

/* CONSTANTES PARA ALGORITMO DE GRAVITACION */
#define PI (3.141592653589793)
#define M_PI (3.14159265358979323846)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 // Hidrogeno molecular

/* ESTRUCTURAS Y VARIABLES PARA ALGORITMO DE GRAVITACION */
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
float *lastPositionX_slice, *lastPositionY_slice, *lastPositionZ_slice;
float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

cuerpo_t *cuerposTotales;
cuerpo_t *cuerpos;
cuerpo_t *tcuerpos;      // Para recibir cuerpos de otros workers
int delta_tiempo = 1.0f; // Intervalo de tiempo, longitud de un paso
int pasos;
int N;

// Fuerzas para cada cuerpo
float *fuerza_totalX_slice;
float *fuerza_totalY_slice;
float *fuerza_totalZ_slice;

int rank;
int T_MPI;
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

/* MAIN */
int main(int argc, char const *argv[])
{
    // MPI
    int idW_MPI;
    MPI_Init(&argc, &argv);                  // Inicializa el ambiente MPI. No debe haber sentencias antes.
    MPI_Comm_rank(MPI_COMM_WORLD, &idW_MPI); // Obtiene el identificador de cada proceso (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &T_MPI);   // Obtiene el número total de procesos (size)

    // ARGUMENTOS
    if (argc < 5)
    {
        printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos> <threads>\n", argv[0]);
        return -1;
    }
    N = atoi(argv[1]);
    dt = atoi(argv[2]);
    pasos = atoi(argv[3]);

    // PTHREADS
    int T_PTHREADS = atoi(argv[4]);
    pthread_t threads[T_PTHREADS];
    pthread_barrier_init(&barrera, NULL, T_PTHREADS);

    for (int i = 0; i < T_PTHREADS; i++)
    {
        thread_ids[i] = i;
        pthread_create(&threads[i], NULL, pfunction, &thread_ids[i]);
    }

    // SIMULACION
    int paso;
    int dest;

    int slice_MPI = N / T_MPI;
    int ini_MPI = idW_MPI * slice_MPI;
    int lim_MPI = ini_MPI + slice_MPI;

    // Reservar memoria para cuerpos de mi bloque de mensaje
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * slice_MPI);
    // Reservar memoria para cuerpos con el tamaño de mi bloque de mensaje
    tcuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * slice_MPI);
    // Reservar memoria para fuerzas con el tamaño de mi bloque de mensaje
    fuerza_totalX_slice = (float *)malloc(sizeof(float) * slice_MPI);
    fuerza_totalY_slice = (float *)malloc(sizeof(float) * slice_MPI);
    fuerza_totalZ_slice = (float *)malloc(sizeof(float) * slice_MPI);
    // Reservar memoria para guardar las ultimas posiciones de los cuerpos
    lastPositionX_slice = (float *)malloc(sizeof(float) * slice_MPI);
    lastPositionY_slice = (float *)malloc(sizeof(float) * slice_MPI);
    lastPositionZ_slice = (float *)malloc(sizeof(float) * slice_MPI);

    // Si soy el proceso master idW_MPI == 0, inicializo los cuerpos
    if (idW_MPI == 0)
    {
        cuerposTotales = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
        fuerza_totalX = (float *)malloc(sizeof(float) * N);
        fuerza_totalY = (float *)malloc(sizeof(float) * N);
        fuerza_totalZ = (float *)malloc(sizeof(float) * N);
        // Allocate memory for last position vectors
        lastPositionX = (float *)malloc(sizeof(float) * N);
        lastPositionY = (float *)malloc(sizeof(float) * N);
        lastPositionZ = (float *)malloc(sizeof(float) * N);

        inicializarCuerpos(cuerposTotales, N);

        // Enviar a cada worker su bloque de cuerpos
        // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
        MPI_Scatter(cuerposTotales, slice_MPI * sizeof(cuerpo_t), MPI_BYTE, cuerpos, slice_MPI * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    }

    for (paso = 0; paso < pasos; paso++)
    {
        // Paso 1: Enviar mis cuerpos a workers con menor idW y calcular fuerzas para mi bloque

        // Envía los cuerpos a los workers con menor idW
        for (dest = 0; dest < T_MPI; dest++)
        {
            if (dest < idW_MPI)
            {
                // Enviar cuerpos a workers con menor idW
                // int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
                MPI_Send(cuerpos, slice_MPI * sizeof(cuerpo_t), MPI_BYTE, dest, 0, MPI_COMM_WORLD);
            }
        }
        // Calcular fuerzas para mi bloque usando pthreads
        pthread_barrier_wait(&barrera);
        calcularFuerzas(ini_MPI, lim_MPI);
        pthread_barrier_wait(&barrera);

        // Paso 2: Recibir cuerpos de workers con mayor idW y calcular fuerzas entre lo recibido y mi bloque.
        // Luego, enviarselas a estos nuevamente.

        for (dest = idW_MPI + 1; dest < T_MPI; dest++)
        {
            // Recibir cuerpos de workers con mayor idW
            // int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
            MPI_Recv(tcuerpos, slice_MPI * sizeof(cuerpo_t), MPI_BYTE, dest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Calcular fuerzas entre los cuerpos recibidos y mi bloque usando pthreads
            pthread_barrier_wait(&barrera);
            calcularFuerzas(ini_MPI, lim_MPI);
            pthread_barrier_wait(&barrera);

            // Enviar fuerzas slice calculadas a los workers con mayor idW
            MPI_Send(fuerza_totalX_slice, slice_MPI * sizeof(float), MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
            MPI_Send(fuerza_totalY_slice, slice_MPI * sizeof(float), MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
            MPI_Send(fuerza_totalZ_slice, slice_MPI * sizeof(float), MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
        }

        // Paso 3: Recibir fuerzas de los workers con menor idW y actualizar p y v
        float *temp_forcesX = (float *)malloc(sizeof(float) * slice_MPI);
        float *temp_forcesY = (float *)malloc(sizeof(float) * slice_MPI);
        float *temp_forcesZ = (float *)malloc(sizeof(float) * slice_MPI);

        for (dest = 0; dest < idW_MPI; dest++)
        {
            // Recibir fuerzas calculadas por workers con menor idW
            MPI_Recv(temp_forcesX, slice_MPI * sizeof(float), MPI_FLOAT, dest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(temp_forcesY, slice_MPI * sizeof(float), MPI_FLOAT, dest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(temp_forcesZ, slice_MPI * sizeof(float), MPI_FLOAT, dest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Sumar las fuerzas recibidas a las fuerzas totales
            for (int i = 0; i < slice_MPI; i++)
            {
                fuerza_totalX_slice[i] += temp_forcesX[i];
                fuerza_totalY_slice[i] += temp_forcesY[i];
                fuerza_totalZ_slice[i] += temp_forcesZ[i];
            }
        }

        // Actualizar posiciones y velocidades usando pthreads
        pthread_barrier_wait(&barrera);
        moverCuerpos(ini_MPI, lim_MPI);
        pthread_barrier_wait(&barrera);

        free(temp_forcesX);
        free(temp_forcesY);
        free(temp_forcesZ);

        // Paso 4: Reinicializar f a cero
        for (int i = 0; i < slice_MPI; i++)
        {
            fuerza_totalX_slice[i] = 0.0f;
            fuerza_totalY_slice[i] = 0.0f;
            fuerza_totalZ_slice[i] = 0.0f;
        }
    }

    // Si soy el proceso master idW_MPI == 0, muestro el tiempo total de la simulación
    if (idW_MPI == 0)
    {
        tFin = dwalltime();
        tTotal = tFin - tIni;
        printf("Tiempo en segundos: %f\n", tTotal);

        // Recibir con Gather las ultimas posiciones de los cuerpos de todos los workers
        // int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
        MPI_Gather(cuerpos, slice_MPI * sizeof(cuerpo_t), MPI_BYTE, cuerposTotales, slice_MPI * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
        // Print last positions of all bodies
        printf("\n=== Last Positions of Bodies ===\n");
        printf("%-6s %-15s %-15s %-15s\n", "ID", "X", "Y", "Z");
        for (int i = 0; i < N; i++)
        {
            printf("%-6d %-15.6f %-15.6f %-15.6f\n", i, cuerposTotales[i].px, cuerposTotales[i].py, cuerposTotales[i].pz);
        }
    }

    finalizar();
    MPI_Finalize();
    return 0;
}

/* FUNCIONES PARA TIEMPOS */
double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

/* FUNCIONES PARA PTHREADS */
void *pfunction(void *arg)
{
    int idW_PTHREADS = *((int *)arg);
    int paso;
    // printf("Hello from thread %d of %d\n", idW, T);
    int slice = N / T;
    int ini = idW_PTHREADS * slice;
    int lim = ini + slice;

    for (paso = 0; paso < pasos; paso++)
    {
        calcularFuerzas(ini, lim);
        pthread_barrier_wait(&barrera);
        moverCuerpos(ini, lim);
        pthread_barrier_wait(&barrera);
    }
    pthread_exit(NULL);
}


/* FUNCIONES DE GRAVITACION */
void calcularFuerzas(int ini, int lim)
{
    int cuerpo1, cuerpo2;
    float dif_X, dif_Y, dif_Z;
    float distancia;
    float F;

    for (cuerpo1 = ini; cuerpo1 < lim; cuerpo1++)
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
            if (cuerpo2 >= ini && cuerpo2 < lim)
            {
                fuerza_totalX[cuerpo2] -= dif_X;
                fuerza_totalY[cuerpo2] -= dif_Y;
                fuerza_totalZ[cuerpo2] -= dif_Z;
            }
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


/* PROCESOS DE INICIALIZACION DE CUERPOS */
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