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

int dt = 1.0f;
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
void calcularFuerzasInternas(cuerpo_t *bloque_cuerpos, double *F1_totalX, double *F1_totalY, double *F1_totalZ);
void calcularFuerzasEntreBloques(cuerpo_t *bloque_cuerpos_1, cuerpo_t *bloque_cuerpos_2, double *F1_totalX, double *F1_totalY, double *F1_totalZ, double *F2_totalX, double *F2_totalY, double *F2_totalZ);
void moverCuerpos(cuerpo_t *cuerpos_totales, double *F_totalX, double *F_totalY, double *F_totalZ, int blocksize);
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
    dt = atof(argv[2]);
    pasos = atoi(argv[3]);
    

    // Lógica del worker para cada proceso
    Worker();

    // Fin de MPI
    MPI_Finalize();
    return 0;
}


void Worker(void)
{
    //Vector para todos los cuerpos
    cuerpo_t *cuerpos_totales = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);

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
        // Inicializamos todos los cuerpos
        inicializarCuerpos(cuerpos_totales, N);


        // Copiamos la primera mitad de cuerpos al bloque local del worker
        memcpy(cuerpos_local, cuerpos_totales, blockSize * sizeof(cuerpo_t));

        // Enviamos la segunda mitad de cuerpos al worker 1
        MPI_Send(&cuerpos_totales[blockSize], blockSize * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD);
    }
    else
    { 
        // Proceso 1 recibe su bloque inicial de cuerpos del Proceso 0 (cero)
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

        // Calculo las fuerzas dentro de mi bloque
        calcularFuerzasInternas(cuerpos_local, fuerzas_localX, fuerzas_localY, fuerzas_localZ);
     

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
            
            // Calcular fuerzas entre mi bloque y el bloque recibido de worker 1

            calcularFuerzasEntreBloques(cuerpos_local, cuerpos_temp, fuerzas_localX, fuerzas_localY, fuerzas_localZ, fuerzas_tempX, fuerzas_tempY, fuerzas_tempZ);

            
            // Enviar las fuerzas calculadas (B0 sobre B1) al Worker 1
            MPI_Send(fuerzas_tempX, tempSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Send(fuerzas_tempY, tempSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Send(fuerzas_tempZ, tempSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }

        // Etapa 3: Recibir las fuerzas de workers con menor idW y Actualizar p y v 
        // Solo Worker 1 recibe del Worker 0 
        if (idW_MPI == 1)
        { 
            
            // Recibe las fuerzas de B0 sobre B1 del Worker 0
            MPI_Recv(fuerzas_tempX, blockSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(fuerzas_tempY, blockSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(fuerzas_tempZ, blockSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            // Combina las fuerzas recibidas con las propias
            for (int i = 0; i < blockSize; i++)
            {
                fuerzas_localX[i] += fuerzas_tempX[i];
                fuerzas_localY[i] += fuerzas_tempY[i];
                fuerzas_localZ[i] += fuerzas_tempZ[i];
            }

        }

        // Cada cuerpo actualiza sus posiciones y velocidades 
        moverCuerpos(cuerpos_local, fuerzas_localX, fuerzas_localY, fuerzas_localZ, blockSize);

        if(idW_MPI == 1){
            // Envio mis cuerpos actualizados para que el worker cero actualice el vector total de cuerpos
            MPI_Send(cuerpos_local, blockSize * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }
        else{
            // Worker cero actualiza sus cuerpos locales dentro del vector total
            memcpy(cuerpos_totales, cuerpos_local, blockSize * sizeof(cuerpo_t));

            // Recibe de worker 1 sus cuerpos actualizados para colocarlos dentro del vector total
            MPI_Recv(&cuerpos_totales[blockSize], blockSize * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Todos los workers reciben el vector total actualizado
        MPI_Bcast(cuerpos_totales, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);

        // Cada worker copia su parte correspondiente de cuerpos locales
        memcpy(cuerpos_local, &cuerpos_totales[idW_MPI * blockSize], blockSize * sizeof(cuerpo_t));
    }


    tFin = dwalltime();
    tTotal = tFin - tIni;


    // ====================================================================
    // RECOLECCIÓN FINAL DE RESULTADOS (SOLO PROCESO 0 IMPRIME)
    // ====================================================================

    if(idW_MPI == 0)
    { 

        printf("Tiempo en segundos: %f\n", tTotal);
        for (int i = 0; i < N; i++)
        {
            printf("%f\n%f\n%f\n", cuerpos_totales[i].px, cuerpos_totales[i].py, cuerpos_totales[i].pz);
        }
    }

    free(cuerpos_totales);
    free(cuerpos_local);
    free(fuerzas_localX);
    free(fuerzas_localY);
    free(fuerzas_localZ);
    free(cuerpos_temp);
    free(fuerzas_tempX);
    free(fuerzas_tempY);
    free(fuerzas_tempZ);
}

void calcularFuerzasInternas(cuerpo_t *bloque_cuerpos, double *F1_totalX, double *F1_totalY, double *F1_totalZ)
{
    int cuerpo1, cuerpo2;
    double dif_X, dif_Y, dif_Z;
    double distancia;
    double F;

    int slice = N / T_MPI;

    for (cuerpo1 = 0; cuerpo1 < slice - 1; cuerpo1++)
    {
        for (cuerpo2 = cuerpo1 + 1; cuerpo2 < slice; cuerpo2++)
        {
            if ((bloque_cuerpos[cuerpo1].px == bloque_cuerpos[cuerpo2].px) &&
                (bloque_cuerpos[cuerpo1].py == bloque_cuerpos[cuerpo2].py) &&
                (bloque_cuerpos[cuerpo1].pz == bloque_cuerpos[cuerpo2].pz))
                continue;

            dif_X = bloque_cuerpos[cuerpo2].px - bloque_cuerpos[cuerpo1].px;
            dif_Y = bloque_cuerpos[cuerpo2].py - bloque_cuerpos[cuerpo1].py;
            dif_Z = bloque_cuerpos[cuerpo2].pz - bloque_cuerpos[cuerpo1].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            F = (G * bloque_cuerpos[cuerpo1].masa * bloque_cuerpos[cuerpo2].masa) / (distancia * distancia);

            F1_totalX[cuerpo1] += dif_X * F;
            F1_totalY[cuerpo1] += dif_Y * F;
            F1_totalZ[cuerpo1] += dif_Z * F;

            F1_totalX[cuerpo2] -= dif_X * F;
            F1_totalY[cuerpo2] -= dif_Y * F;
            F1_totalZ[cuerpo2] -= dif_Z * F;
        }
    }
}

void calcularFuerzasEntreBloques(cuerpo_t *bloque_cuerpos_1, cuerpo_t *bloque_cuerpos_2, double *F1_totalX, double *F1_totalY, double *F1_totalZ, double *F2_totalX, double *F2_totalY, double *F2_totalZ)
{
    int cuerpo1, cuerpo2;
    double dif_X, dif_Y, dif_Z;
    double distancia;
    double F;

    int slice = N / T_MPI;

    for (cuerpo1 = 0; cuerpo1 < slice; cuerpo1++)
    {
        for (cuerpo2 = 0; cuerpo2 < slice; cuerpo2++)
        {
            if ((bloque_cuerpos_1[cuerpo1].px == bloque_cuerpos_2[cuerpo2].px) &&
                (bloque_cuerpos_1[cuerpo1].py == bloque_cuerpos_2[cuerpo2].py) &&
                (bloque_cuerpos_1[cuerpo1].pz == bloque_cuerpos_2[cuerpo2].pz))
                continue;

            dif_X = bloque_cuerpos_2[cuerpo2].px - bloque_cuerpos_1[cuerpo1].px;
            dif_Y = bloque_cuerpos_2[cuerpo2].py - bloque_cuerpos_1[cuerpo1].py;
            dif_Z = bloque_cuerpos_2[cuerpo2].pz - bloque_cuerpos_1[cuerpo1].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            F = (G * bloque_cuerpos_1[cuerpo1].masa * bloque_cuerpos_2[cuerpo2].masa) / (distancia * distancia);

            F1_totalX[cuerpo1] += dif_X * F;
            F1_totalY[cuerpo1] += dif_Y * F;
            F1_totalZ[cuerpo1] += dif_Z * F;

            F2_totalX[cuerpo2] -= dif_X * F;
            F2_totalY[cuerpo2] -= dif_Y * F;
            F2_totalZ[cuerpo2] -= dif_Z * F;
        }
    }
}


void moverCuerpos(cuerpo_t *cuerpos_totales, double *F_totalX, double *F_totalY, double *F_totalZ, int blocksize)
{
    int cuerpo;

    for (cuerpo = 0; cuerpo < blocksize; cuerpo++)
    {

		F_totalX[cuerpo] *= 1 / cuerpos_totales[cuerpo].masa;
		F_totalY[cuerpo] *= 1 / cuerpos_totales[cuerpo].masa;
		// F_totalZ[cuerpo] *= 1/cuerpos_totales[cuerpo].masa;

		cuerpos_totales[cuerpo].vx += F_totalX[cuerpo] * dt;
		cuerpos_totales[cuerpo].vy += F_totalY[cuerpo] * dt;
		// cuerpos_totales[cuerpo].vz += F_totalZ[cuerpo]*dt;

		cuerpos_totales[cuerpo].px += cuerpos_totales[cuerpo].vx * dt;
		cuerpos_totales[cuerpo].py += cuerpos_totales[cuerpo].vy * dt;
		// cuerpos_totales[cuerpo].pz += cuerpos_totales[cuerpo].vz *dt;

        F_totalX[cuerpo] = 0.0;
		F_totalY[cuerpo] = 0.0;
		F_totalZ[cuerpo] = 0.0;
    }



    for (int i = 0; i < N; i++)
    {
        printf("%f\n%f\n%f\n", cuerpos_totales[i].vx, cuerpos_totales[i].vy);
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

