#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

double toroide_alfa;
double toroide_theta;
double toroide_incremento;
double toroide_lado;
double toroide_r;
double toroide_R;

double dt = 1.0f;
int pasos;
int N;

// Variables MPI
int idW_MPI; // ID del proceso MPI
int T_MPI; // Total de procesos (asumimos dos para nuestro caso)

void inicializarEstrella(cuerpo_t *cuerpo, int i, double n);
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n);
void inicializarH2(cuerpo_t *cuerpo, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos_arr, int N_total);

// Funciones para MPI
void calcularFuerzas(cuerpo_t *cuerpo_influencia, cuerpo_t *cuerpo_afectado, double *F_totalX_afectado, double *F_totalY_afectado, double *F_totalZ_afectado);
void moverCuerpos(cuerpo_t *cuerpos_arr, double *fuerzasX_arr, double *fuerzasY_arr, double *fuerzasZ_arr, double dt_val, int num_cuerpos_local);void Worker(void);
void Worker(void);




int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &idW_MPI);
    MPI_Comm_size(MPI_COMM_WORLD, &T_MPI);

    if (argc < 4)
    {
        fprintf(stderr, "Ejecutar: %s <nro. de cuerpos> <DT> <pasos>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    N = atoi(argv[1]);
    dt = atoi(argv[2]);
    pasos = atoi(argv[3]);

    // Lógica del worker para cada proceso
    Worker();

    // Fin de MPI
    MPI_Finalize();
    return 0;
}


void Worker(void)
{
    // Tamaño de mi bloque
    int blockSize = N/T_MPI;

    // Variables locales para cada worker
    cuerpo_t *cuerpos_local; // Mi bloque de cuerpos
    double *fuerzas_localX, *fuerzas_localY, *fuerzas_localZ; // Fuerzas acumuladas para mis cuerpos

    // Buffers temporales para comunicación con otros workers
    int tempSize = blockSize; // Son iguales por asumir dos procesos solamente
    cuerpo_t *cuerpos_temp;                // Buffer para cuerpos recibidos
    double *fuerzas_tempX, *fuerzas_tempY, *fuerzas_tempZ; // Buffer para fuerzas temporales (enviar/recibir)

    // Asignar memoria para las variables locales y buffers
    cuerpos_local = (cuerpo_t *)malloc(sizeof(cuerpo_t) * blockSize);
    fuerzas_localX = (double *)malloc(sizeof(double) * blockSize);
    fuerzas_localY = (double *)malloc(sizeof(double) * blockSize);
    fuerzas_localZ = (double *)malloc(sizeof(double) * blockSize);

    cuerpos_temp = (cuerpo_t *)malloc(sizeof(cuerpo_t) * tempSize);
    fuerzas_tempX = (double *)malloc(sizeof(double) * tempSize);
    fuerzas_tempY = (double *)malloc(sizeof(double) * tempSize);
    fuerzas_tempZ = (double *)malloc(sizeof(double) * tempSize);


    // ====================================================================
    // INICIALIZACION GLOBAL Y DISTRIBUCION INICIAL DE CUERPOS
    // ====================================================================

    // El proceso 0 (cero) realiza la inicialización
    if (idW_MPI == 0)
    {
        cuerpo_t *all_cuerpos_init = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);

        inicializarCuerpos(all_cuerpos_init, N);

        // Copiamos la primera mitad de cuerpos al bloque local del worker
        memcpy(cuerpos_local, all_cuerpos_init, blockSize * sizeof(cuerpo_t));

        // Enviamos la segunda mitad de cuerpos al worker 1
        MPI_Send(&all_cuerpos_init[blockSize], blockSize * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD);

        free(all_cuerpos_init); // Liberamos el espacio para el buffer temporal global
    }
    else
    { // Proceso 1 recibe su bloque inicial de cuerpos del Proceso 0 (cero)
        MPI_Recv(cuerpos_local, blockSize * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Sincronizamos todos los procesos antes de iniciar la simulacion
    MPI_Barrier(MPI_COMM_WORLD);


    // ====================================================================
    // BUCLE PRINCIPAL DE SIMULACION (PASOS DE TIEMPO)
    // ====================================================================

    tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++)
    {
        // Reinicializar f a cero (para mis cuerpos)
        for (int i = 0; i < blockSize; i++)
        {
            fuerzas_localX[i] = 0.0;
            fuerzas_localY[i] = 0.0;
            fuerzas_localZ[i] = 0.0;
        }

        // Etapa 1: Enviar mis cuerpos a workers con menor idW 
        // En este caso solamente proceso 1 es el que debe enviar sus cuerpos
        
        if (idW_MPI == 1)
        {
            MPI_Send(cuerpos_local, blockSize * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }

        // Calcular fuerzas para mi propio bloque (B_0, B_0) 
        for (int i = 0; i < blockSize; i++)
        { // cuerpo1 en mi bloque
            for (int j = 0; j < blockSize; j++)
            { // cuerpo2 en mi bloque
                if (i == j)
                    continue;
                // Fuerza que cuerpos_local[j] ejerce sobre cuerpos_local[i]
                calcularFuerzas(&cuerpos_local[j], &cuerpos_local[i], &fuerzas_localX[i], &fuerzas_localY[i], &fuerzas_localZ[i]);
            }
        }

        // Etapa 2: Recibir cuerpos de workers con mayor idW, calcular fuerzas y enviarlas 
        // Solo Worker 0 (cero) recibe del Worker 1
        if (idW_MPI == 0)
        { 
            // Recibimos cuerpos del Worker 1
            MPI_Recv(cuerpos_temp, tempSize * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Reinicializar fuerzas_temp 
            for (int i = 0; i < tempSize; i++)
            {
                fuerzas_tempX[i] = 0.0;
                fuerzas_tempY[i] = 0.0;
                fuerzas_tempZ[i] = 0.0;
            }

            // Calcular fuerzas entre MI bloque (B0) y el bloque recibido (B1)
            for (int i = 0; i < blockSize; i++)
            { // cuerpo del bloque 0 (B0)
                for (int j = 0; j < tempSize; j++)
                { // cuerpo del bloque 1 (B1)
                    // Fuerza que cuerpos_temp[j] (B1) ejerce sobre cuerpos_local[i] (B0)
                    calcularFuerzas(&cuerpos_temp[j], &cuerpos_local[i], &fuerzas_localX[i], &fuerzas_localY[i], &fuerzas_localZ[i]);
                    // Fuerza que cuerpos_local[i] (B0) ejerce sobre cuerpos_temp[j] (B1)
                    // Esta fuerza es la que le pertenece a Worker 1, por eso se acumula en fuerzas_temp
                    calcularFuerzas(&cuerpos_local[i], &cuerpos_temp[j], &fuerzas_tempX[j], &fuerzas_tempY[j], &fuerzas_tempZ[j]);
                }
            }
            // Enviar las fuerzas calculadas (B0 sobre B1) al Worker 1
            MPI_Send(fuerzas_tempX, tempSize, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            MPI_Send(fuerzas_tempY, tempSize, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
            MPI_Send(fuerzas_tempZ, tempSize, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD);
        }

        // Etapa 3: Recibir las fuerzas de workers con menor idW y Actualizar p y v 
        // Solo Worker 1 recibe del Worker 0 
        if (idW_MPI == 1)
        { 
            // Recibe las fuerzas de B0 sobre B1 del Worker 0
            MPI_Recv(fuerzas_tempX, blockSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(fuerzas_tempY, blockSize, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(fuerzas_tempZ, blockSize, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Combina las fuerzas recibidas con las suyas propias (B1 sobre B1)
            for (int i = 0; i < blockSize; i++)
            {
                fuerzas_localX[i] += fuerzas_tempX[i];
                fuerzas_localY[i] += fuerzas_tempY[i];
                fuerzas_localZ[i] += fuerzas_tempZ[i];
            }
        }

        // Actualizar posiciones y velocidades para mis cuerpos (cuerpos_local)
        moverCuerpos(cuerpos_local, fuerzas_localX, fuerzas_localY, fuerzas_localZ, dt, blockSize);
    }

    tFin = dwalltime();
    tTotal = tFin - tIni;


    // ====================================================================
    // RECOLECCIÓN FINAL DE RESULTADOS (SOLO PROCESO 0 IMPRIME)
    // ====================================================================

    if (idW_MPI == 0)
    { 
        cuerpo_t *final_bodies = (cuerpo_t *)malloc(N * sizeof(cuerpo_t));

        // Copia sus propios cuerpos finales a la primera mitad de final_bodies
        memcpy(final_bodies, cuerpos_local, blockSize * sizeof(cuerpo_t));

        // Recibe los cuerpos finales del Worker 1 y los coloca en la segunda mitad
        MPI_Recv(&final_bodies[blockSize], blockSize * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Tiempo en segundos: %f\n", tTotal);
        for (int i = 0; i < N; i++)
        {
            printf("%f\n%f\n%f\n", final_bodies[i].px, final_bodies[i].py, final_bodies[i].pz);
        }
        free(final_bodies);
    }
    else
    { // Worker 1 envía sus cuerpos finales a Worker 0
        MPI_Send(cuerpos_local, blockSize * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }

    free(cuerpos_local);
    free(fuerzas_localX);
    free(fuerzas_localY);
    free(fuerzas_localZ);
    free(cuerpos_temp);
    free(fuerzas_tempX);
    free(fuerzas_tempY);
    free(fuerzas_tempZ);
}

// Calcula la fuerza que 'cuerpo_influencia' ejerce sobre 'cuerpo_afectado'
// y la acumula en los punteros de fuerza total de 'cuerpo_afectado'.
void calcularFuerzas(cuerpo_t *cuerpo_influencia, cuerpo_t *cuerpo_afectado, double *F_totalX_afectado, double *F_totalY_afectado, double *F_totalZ_afectado)
{
    if ((cuerpo_influencia->px == cuerpo_afectado->px) && (cuerpo_influencia->py == cuerpo_afectado->py) && (cuerpo_influencia->pz == cuerpo_afectado->pz))
        return;

    double dif_X = cuerpo_afectado->px - cuerpo_influencia->px;
    double dif_Y = cuerpo_afectado->py - cuerpo_influencia->py;
    double dif_Z = cuerpo_afectado->pz - cuerpo_influencia->pz;

    double distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

    double F_magnitud = (G * cuerpo_influencia->masa * cuerpo_afectado->masa) / (distancia * distancia);

    // Componentes de la fuerza (Fx, Fy, Fz)
    *F_totalX_afectado += (dif_X / distancia) * F_magnitud;
    *F_totalY_afectado += (dif_Y / distancia) * F_magnitud;
    *F_totalZ_afectado += (dif_Z / distancia) * F_magnitud;
}

// Mueve los cuerpos en el arreglo 'cuerpos_arr' usando las fuerzas acumuladas
// en 'fuerzasX_arr/Y_arr/Z_arr'. Opera sobre 'num_cuerpos_local'.
void moverCuerpos(cuerpo_t *cuerpos_arr, double *fuerzasX_arr, double *fuerzasY_arr, double *fuerzasZ_arr, double dt_val, int num_cuerpos_local)
{
    for (int i = 0; i < num_cuerpos_local; i++)
    {
        // Aceleracion = Fuerza / Masa
        double ax = fuerzasX_arr[i] / cuerpos_arr[i].masa;
        double ay = fuerzasY_arr[i] / cuerpos_arr[i].masa;
    //  double az = fuerzasZ_arr[i] / cuerpos_arr[i].masa;

        // Actualizar velocidad
        cuerpos_arr[i].vx += ax * dt_val;
        cuerpos_arr[i].vy += ay * dt_val;
    //  cuerpos_arr[i].vz += az * dt_val;

        // Actualizar posicion
        cuerpos_arr[i].px += cuerpos_arr[i].vx * dt_val;
        cuerpos_arr[i].py += cuerpos_arr[i].vy * dt_val;
    //  cuerpos_arr[i].pz += cuerpos_arr[i].vz * dt_val;
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

