#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <mpi.h>

/* Para tiempo de ejecucion */
double dwalltime();
double tIni, tFin, tTotal;

/* Constantes para Algoritmo de gravitacion */
#define PI (3.141592653589793)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 // Hidrogeno molecular

/* Pthreads */
void *pfunction(void *arg);
pthread_barrier_t barrier_threads; // Barrera global
pthread_barrier_t barrier_main;    // Barrera para el main

/* Argumentos */
int N;          // Número de cuerpos
int dt;         // Intervalo de tiempo, longitud de un paso
int pasos;      // Número de pasos a simular
int T_MPI;      // Número total de procesos MPI
int T_PTHREADS; // Número de threads por proceso
int idW_MPI;    // Identificador del proceso MPI
/* Estructuras y variables para Algoritmo de gravitacion */
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
float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

/* Simulacion */
void calcularFuerzas(int ini, int lim, int lim_block);
void calcularFuerzasEntreBloques(int bloque1_ini, int bloque1_fin, int bloque2_ini, int bloque2_fin);
void moverCuerpos(int ini, int lim);
cuerpo_t *cuerpos;
cuerpo_t *recv_cuerpos;
float *recv_fuerza_totalX, *recv_fuerza_totalY, *recv_fuerza_totalZ;
float *recv_fuerza_totalXYZ;

/* Inicializacion de Cuerpos */
void inicializarEstrella(cuerpo_t *cuerpo, int i, double n);
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n);
void inicializarH2(cuerpo_t *cuerpo, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos, int N);

/* Finalizar */
void finalizar(void);

/* Variables globales para threads */
int bloque2_ini_global;  // Para calcularFuerzasEntreBloqueThread
int bloque2_fin_global;  // Para calcularFuerzasEntreBloqueThread
void *calcularFuerzasThread(void *arg);
void *calcularFuerzasEntreBloqueThread(void *arg);
void *moverCuerposThread(void *arg);

/* MAIN */
int main(int argc, char *argv[])
{
    // MPI
    int mpi_error = MPI_Init(&argc, &argv); // Check MPI initialization
    if (mpi_error != MPI_SUCCESS)
    {
        fprintf(stderr, "Error initializing MPI.\n");
        return -1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &idW_MPI); // Obtiene el identificador de cada proceso (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &T_MPI);   // Obtiene el número total de procesos (size)

    // ARGUMENTOS
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

    if (N <= 0 || dt <= 0 || pasos <= 0 || T_PTHREADS <= 0)
    {
        fprintf(stderr, "Error: All arguments must be positive integers.\n");
        MPI_Finalize();
        return -1;
    }



    int slice_MPI = N / T_MPI;         // Número de cuerpos por proceso
    int ini_MPI = idW_MPI * slice_MPI; // Índice inicial para este proceso
    int lim_MPI = ini_MPI + slice_MPI; // Índice final para este proceso

    // Inicializar el arreglo global
    for (int i = 0; i < 512; i++)
    {
        array[i] = i;
        recv_array[i] = 0;
    }

    // Inicializar cuerpos y fuerzas
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (float *)malloc(sizeof(float) * N);
    fuerza_totalY = (float *)malloc(sizeof(float) * N);
    fuerza_totalZ = (float *)malloc(sizeof(float) * N);
    recv_cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    recv_fuerza_totalX = (float *)malloc(sizeof(float) * N);
    recv_fuerza_totalY = (float *)malloc(sizeof(float) * N);
    recv_fuerza_totalZ = (float *)malloc(sizeof(float) * N);

    recv_fuerza_totalXYZ = (float *)malloc(sizeof(float) * N * 3); // Para recibir fuerzas en X, Y, Z

    if (!cuerpos || !fuerza_totalX || !fuerza_totalY || !fuerza_totalZ || !recv_cuerpos || !recv_fuerza_totalX || !recv_fuerza_totalY || !recv_fuerza_totalZ)
    {
        fprintf(stderr, "Error allocating memory.\n");
        MPI_Finalize();
        return -1;
    }
    // Inicializar cuerpos
    inicializarCuerpos(cuerpos, N);

    MPI_Barrier(MPI_COMM_WORLD);



    // Inicializar la barrera para T_PTHREADS hilos
    pthread_barrier_init(&barrier_threads, NULL, T_PTHREADS);

    if (idW_MPI == 0)
    {

        tIni = dwalltime(); // Iniciar tiempo de simulación
    }

    for (int paso = 0; paso < pasos; paso++)
    {
        /* ====== */
        /* PASO 1 */
        /* ====== */



        // Enviar mi arreglo de cuerpos a todos los procesos MPI con menor id
        for (int i = 0; i < T_MPI; i++)
        {
            if (i < idW_MPI)
            {
                // Calculate exact size to send
                MPI_Status status;
                int send_size = slice_MPI * sizeof(cuerpo_t);

                if (MPI_Send(cuerpos + ini_MPI, send_size, MPI_BYTE, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
                {
                    fprintf(stderr, "Error enviando cuerpos de proceso %d a proceso %d\n", idW_MPI, i);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }



        // Crear y ejecutar threads para calcularFuerzas
        pthread_t threads[T_PTHREADS];
        int thread_ids[T_PTHREADS];
        
        for(int i = 0; i < T_PTHREADS; i++) {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, calcularFuerzasThread, (void *)&thread_ids[i]);
        }
        
        for(int i = 0; i < T_PTHREADS; i++) {
            pthread_join(threads[i], NULL);
        }

        /* ====== */
        /* PASO 2 */
        /* ====== */



        // Recibir los cuerpos de procesos MPI con mayor id en recv_cuerpos
        for (int i = idW_MPI + 1; i < T_MPI; i++)
        {
            MPI_Status status;
            // Calculate exact size to receive
            int recv_size = slice_MPI * sizeof(cuerpo_t);

            if (MPI_Recv(recv_cuerpos + (i * slice_MPI), recv_size, MPI_BYTE, i, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
            {
                fprintf(stderr, "Error recibiendo cuerpos de proceso %d en proceso %d\n", i, idW_MPI);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // Verify received size
            int received_count;
            MPI_Get_count(&status, MPI_BYTE, &received_count);
            if (received_count != recv_size)
            {
                fprintf(stderr, "Error: mensaje truncado. Esperados %d bytes, recibidos %d\n",
                        recv_size, received_count);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }



            // Actualizar solo el slice recibido en el arreglo cuerpos
            for (int j = 0; j < slice_MPI; j++)
            {
                cuerpos[j + (i * slice_MPI)] = recv_cuerpos[j + (i * slice_MPI)];
            }
            
            // Establecer variables globales para calcularFuerzasEntreBloques
            bloque2_ini_global = i * slice_MPI;
            bloque2_fin_global = N;
            
            // Crear y ejecutar threads para calcularFuerzasEntreBloques
            for(int i = 0; i < T_PTHREADS; i++) {
                thread_ids[i] = i;
                pthread_create(&threads[i], NULL, calcularFuerzasEntreBloqueThread, (void *)&thread_ids[i]);
            }
            
            for(int i = 0; i < T_PTHREADS; i++) {
                pthread_join(threads[i], NULL);
            }

            /*
            FALTA: "calculo fuerzas entre mi bloque y el bloque de otherWorker, los guardo en tf"
            FALTA: "send forces[otherWorker](tf)"
            */

            // Guardo fuerzas X Y Z en XYZ
            for (int x = 0; x < N; x++)
            {
                recv_fuerza_totalXYZ[x * 3] = fuerza_totalX[x];
                recv_fuerza_totalXYZ[x * 3 + 1] = fuerza_totalY[x];
                recv_fuerza_totalXYZ[x * 3 + 2] = fuerza_totalZ[x];
            }

            // Enviar fuerzas total XYZ calculadas a procesos MPI con MAYOR id
            if (i > idW_MPI)
            {
                // Send exactly 3 floats per body in the slice
                int send_size = slice_MPI * 3; // Number of float values to send

                // Send the forces for this slice
                if (MPI_Send(recv_fuerza_totalXYZ + (ini_MPI * 3), send_size, MPI_FLOAT, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
                {
                    fprintf(stderr, "Error sending forces XYZ from process %d to process %d\n", idW_MPI, i);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }


        }

        /* ====== */
        /* PASO 3 */
        /* ====== */



        // Recibir fuerzas de procesos MPI con menor id
        for (int i = 0; i < idW_MPI; i++)
        {
            MPI_Status status;
            int recv_size = slice_MPI * 3; // Number of float values to receive

            // Receive forces for this slice
            if (MPI_Recv(recv_fuerza_totalXYZ + (i * slice_MPI * 3), recv_size, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
            {
                fprintf(stderr, "Error receiving forces XYZ from process %d in process %d\n", i, idW_MPI);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // Verify received size
            int received_count;
            MPI_Get_count(&status, MPI_FLOAT, &received_count);
            if (received_count != recv_size)
            {
                fprintf(stderr, "Error: truncated message. Expected %d floats, received %d\n",
                        recv_size, received_count);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // Add received forces to total forces arrays
            for (int j = 0; j < slice_MPI; j++)
            {
                int idx = j + (i * slice_MPI);
                fuerza_totalX[idx] += recv_fuerza_totalXYZ[(idx * 3)];
                fuerza_totalY[idx] += recv_fuerza_totalXYZ[(idx * 3) + 1];
                fuerza_totalZ[idx] += recv_fuerza_totalXYZ[(idx * 3) + 2];
            }


        }
        // MOVEMOS LOS CUERPOS

        
        // Crear y ejecutar threads para moverCuerpos
        for(int i = 0; i < T_PTHREADS; i++) {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, moverCuerposThread, (void *)&thread_ids[i]);
        }
        
        for(int i = 0; i < T_PTHREADS; i++) {
            pthread_join(threads[i], NULL);
        }

        /* ====== */
        /* PASO 4 */
        /* ====== */



        // Reinicizar fuerzas totales
        for (int i = 0; i < N; i++)
        {
            fuerza_totalX[i] = 0.0;
            fuerza_totalY[i] = 0.0;
            fuerza_totalZ[i] = 0.0;
        }
        // Reiniciar fuerzas_totalXYZ
        for (int i = 0; i < N * 3; i++)
        {
            recv_fuerza_totalXYZ[i] = 0.0;
        }
    }

    
    // Destruir la barrera
    pthread_barrier_destroy(&barrier_threads);

    MPI_Finalize();

    // Print last positions of all bodies using cuerpos
    if (idW_MPI == 0)
    {
        printf("\n=== Last Positions of Bodies ===\n");
        printf("%-6s %-15s %-15s %-15s\n", "ID", "X", "Y", "Z");
        for (int i = 0; i < N; i++)
        {
            printf("%-6d %-15.6f %-15.6f %-15.6f\n", i, cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
        }
    }

    printf("\n[Proceso %d] ========== SIMULACIÓN COMPLETADA ==========\n", idW_MPI);
    if (idW_MPI == 0)
    {
        tFin = dwalltime(); // Finalizar tiempo de simulación
        tTotal = tFin - tIni;
        printf("Tiempo total de simulación: %f segundos\n", tTotal);
    }

    return 0;
}

/* Para tiempo de ejecucion */
double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

/* Pthreads */

/* Simulacion */ /*  inicio -  final*/
void calcularFuerzas(int ini, int lim, int lim_block)
{
    int cuerpo1, cuerpo2;
    float dif_X, dif_Y, dif_Z;
    float distancia;
    float F;

    for (cuerpo1 = ini; cuerpo1 < lim; cuerpo1++)
    {
        for (cuerpo2 = cuerpo1 + 1; cuerpo2 < lim_block; cuerpo2++)
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

/* Inicializacion de Cuerpos */
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

/* Finalizar */
void finalizar(void)
{
    free(cuerpos);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
    free(recv_cuerpos);
    free(recv_fuerza_totalX);
    free(recv_fuerza_totalY);
    free(recv_fuerza_totalZ);
}

// Añadir esta nueva función:
void calcularFuerzasEntreBloques(int bloque1_ini, int bloque1_fin, int bloque2_ini, int bloque2_fin)
{
    float dif_X, dif_Y, dif_Z;
    float distancia;
    float F;

    // Calcular fuerzas solo entre los cuerpos del bloque1 y bloque2
    for (int cuerpo1 = bloque1_ini; cuerpo1 < bloque1_fin; cuerpo1++)
    {
        for (int cuerpo2 = bloque2_ini; cuerpo2 < bloque2_fin; cuerpo2++)
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

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            // Actualizar fuerzas para el cuerpo1 (siempre en mi bloque local)
            fuerza_totalX[cuerpo1] += dif_X;
            fuerza_totalY[cuerpo1] += dif_Y;
            fuerza_totalZ[cuerpo1] += dif_Z;

            // Actualizar fuerzas para el cuerpo2 (del bloque remoto)
            // Estas fuerzas serán enviadas de vuelta al proceso propietario
            recv_fuerza_totalX[cuerpo2] -= dif_X;
            recv_fuerza_totalY[cuerpo2] -= dif_Y;
            recv_fuerza_totalZ[cuerpo2] -= dif_Z;
        }
    }
}

/* Thread Functions */
void *calcularFuerzasThread(void *arg) {
    int idW = *((int *)arg);
    int slice_MPI = N / T_MPI;
    int ini_MPI = idW_MPI * slice_MPI;
    
    int slice_thread = slice_MPI / T_PTHREADS;
    int ini_thread = ini_MPI + idW * slice_thread;
    int lim_thread = ini_thread + slice_thread;
    
    calcularFuerzas(ini_thread, lim_thread, slice_MPI);
    
    pthread_barrier_wait(&barrier_threads);
    return NULL;
}

void *calcularFuerzasEntreBloqueThread(void *arg) {
    int idW = *((int *)arg);
    int slice_MPI = N / T_MPI;
    int ini_MPI = idW_MPI * slice_MPI;
    
    int slice_thread = slice_MPI / T_PTHREADS;
    int ini_thread = ini_MPI + idW * slice_thread;
    int lim_thread = ini_thread + slice_thread;
    
    // bloque2_ini y bloque2_fin son variables globales que deben ser establecidas antes de crear los threads
    calcularFuerzasEntreBloques(ini_thread, lim_thread, bloque2_ini_global, bloque2_fin_global);
    
    pthread_barrier_wait(&barrier_threads);
    return NULL;
}

void *moverCuerposThread(void *arg) {
    int idW = *((int *)arg);
    int slice_MPI = N / T_MPI;
    int ini_MPI = idW_MPI * slice_MPI;
    
    int slice_thread = slice_MPI / T_PTHREADS;
    int ini_thread = ini_MPI + idW * slice_thread;
    int lim_thread = ini_thread + slice_thread;
    
    moverCuerpos(ini_thread, lim_thread);
    
    pthread_barrier_wait(&barrier_threads);
    return NULL;
}