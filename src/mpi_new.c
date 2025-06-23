#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

double dwalltime();
double tIni, tFin, tTotal;

#define PI (3.141592653589793)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2

// Enum for MPI message tags
enum tag
{
    CUERPOS,
    FUERZAS_X,
    FUERZAS_Y,
    FUERZAS_Z
};

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

// Global variables for force accumulation
double *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;

// Global variables for body initialization
double toroide_alfa;
double toroide_theta;
double toroide_incremento;
double toroide_lado;
double toroide_r;
double toroide_R;

cuerpo_t *cuerpos; // All bodies (on Coordinator)
double dt = 1.0f;
int pasos;
int N; // Total number of bodies

// MPI specific global variables
int idW_MPI; // MPI rank
int T_MPI;   // Total MPI processes
int block_size; // Number of bodies handled by each process
int start_index; // Starting index for bodies for the current process

// Force buffers for MPI communication
double *fuerza_externaX, *fuerza_externaY, *fuerza_externaZ;

// Function declarations
void inicializarEstrella(cuerpo_t *cuerpo, int i, double n);
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n);
void inicializarH2(cuerpo_t *cuerpo, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos, int N);
void finalizar(void);

void calcularFuerzas(cuerpo_t *local_cuerpos, int local_N, cuerpo_t *all_cuerpos, int total_N, int offset);
void moverCuerpos(cuerpo_t *local_cuerpos, int local_N, double delta_tiempo);

void Coordinator(void);
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
    dt = atof(argv[2]); // Using atof for double precision
    pasos = atoi(argv[3]);

    // Distribute bodies (approximately)
    block_size = N / T_MPI;
    // For simplicity, let's make block_size consistent for all processes,
    // and the last process handles any remainder.
    if (idW_MPI == T_MPI - 1) {
        block_size = N - (idW_MPI * (N / T_MPI));
    }

    start_index = idW_MPI * (N / T_MPI);


    // Allocate memory for global force arrays
    // These will be used to accumulate forces for all bodies on each process
    // and will be sent/received.
    fuerza_totalX = malloc(sizeof(double) * N);
    fuerza_totalY = malloc(sizeof(double) * N);
    fuerza_totalZ = malloc(sizeof(double) * N);

    // Allocate memory for external forces (forces from other processes)
    fuerza_externaX = malloc(sizeof(double) * N);
    fuerza_externaY = malloc(sizeof(double) * N);
    fuerza_externaZ = malloc(sizeof(double) * N);

    if (idW_MPI == 0)
    {
        Coordinator();
    }
    else
    {
        Worker();
    }

    // Free global force buffers
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
    free(fuerza_externaX);
    free(fuerza_externaY);
    free(fuerza_externaZ);


    MPI_Finalize();
    return 0;
}

void Coordinator(void)
{
    // Coordinator has all bodies initially
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    inicializarCuerpos(cuerpos, N);

    tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++)
    {
        // 1. Broadcast all body positions to all workers
        MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);

        // 2. Clear force accumulators
        memset(fuerza_totalX, 0, sizeof(double) * N);
        memset(fuerza_totalY, 0, sizeof(double) * N);
        memset(fuerza_totalZ, 0, sizeof(double) * N);

        // 3. Calculate forces for its own assigned block of bodies
        // The Coordinator calculates forces for bodies from index 0 up to block_size-1
        calcularFuerzas(cuerpos, block_size, cuerpos, N, 0);

        // 4. Receive forces from all other workers
        for (int source_rank = 1; source_rank < T_MPI; source_rank++)
        {
            int current_block_size = (source_rank == T_MPI - 1) ? N - (source_rank * (N / T_MPI)) : (N / T_MPI);
            int current_start_index = source_rank * (N / T_MPI);

            MPI_Recv(&fuerza_externaX[current_start_index], current_block_size, MPI_DOUBLE, source_rank, FUERZAS_X, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&fuerza_externaY[current_start_index], current_block_size, MPI_DOUBLE, source_rank, FUERZAS_Y, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&fuerza_externaZ[current_start_index], current_block_size, MPI_DOUBLE, source_rank, FUERZAS_Z, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Accumulate forces from workers into the total force arrays
            for (int i = 0; i < current_block_size; i++) {
                fuerza_totalX[current_start_index + i] += fuerza_externaX[current_start_index + i];
                fuerza_totalY[current_start_index + i] += fuerza_externaY[current_start_index + i];
                fuerza_totalZ[current_start_index + i] += fuerza_externaZ[current_start_index + i];
            }
        }

        // 5. Move all bodies based on the accumulated total forces
        moverCuerpos(cuerpos, N, dt);
    }

    tFin = dwalltime();
    tTotal = tFin - tIni;
    printf("Tiempo en segundos: %f\n", tTotal);

    // Print final positions
    for (int i = 0; i < N; i++)
    {
        printf("%f\n%f\n%f\n", cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
    }

    finalizar(); // Frees 'cuerpos' (all bodies) on Coordinator
}

void Worker(void)
{
    // Workers need a local copy of all bodies to calculate interactions
    cuerpo_t *all_cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    
    // Workers only compute forces for their assigned block of bodies.
    // However, they need to read all body positions to calculate interactions
    // with bodies outside their block.
    cuerpo_t *local_cuerpos_segment = (cuerpo_t *)malloc(sizeof(cuerpo_t) * block_size);

    for (int paso = 0; paso < pasos; paso++)
    {
        // 1. Receive all body positions from the Coordinator
        MPI_Bcast(all_cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);

        // Copy the relevant segment of bodies for this worker
        memcpy(local_cuerpos_segment, all_cuerpos + start_index, sizeof(cuerpo_t) * block_size);

        // 2. Clear force accumulators for the local segment
        // We only care about the forces on our local segment, so we clear only those relevant parts
        // of the global force arrays.
        memset(&fuerza_totalX[start_index], 0, sizeof(double) * block_size);
        memset(&fuerza_totalY[start_index], 0, sizeof(double) * block_size);
        memset(&fuerza_totalZ[start_index], 0, sizeof(double) * block_size);

        // 3. Calculate forces for its own assigned block of bodies
        // local_cuerpos_segment contains the bodies this worker is responsible for.
        // all_cuerpos contains all bodies for interaction calculations.
        // The offset is `start_index` because the force_totalX/Y/Z arrays are for all N bodies,
        // and this worker is responsible for forces from `start_index` to `start_index + block_size - 1`.
        calcularFuerzas(local_cuerpos_segment, block_size, all_cuerpos, N, start_index);

        // 4. Send its calculated forces for its block back to the Coordinator
        MPI_Send(&fuerza_totalX[start_index], block_size, MPI_DOUBLE, 0, FUERZAS_X, MPI_COMM_WORLD);
        MPI_Send(&fuerza_totalY[start_index], block_size, MPI_DOUBLE, 0, FUERZAS_Y, MPI_COMM_WORLD);
        MPI_Send(&fuerza_totalZ[start_index], block_size, MPI_DOUBLE, 0, FUERZAS_Z, MPI_COMM_WORLD);
    }
    
    free(all_cuerpos);
    free(local_cuerpos_segment);
}


// Modified calculateForces to take the local bodies and the full set of bodies
void calcularFuerzas(cuerpo_t *local_cuerpos, int local_N, cuerpo_t *all_cuerpos, int total_N, int offset)
{
    int cuerpo1_idx, cuerpo2_idx; // Indices for the *original* array (all_cuerpos)
    double dif_X, dif_Y, dif_Z;
    double distancia;
    double F;

    // This loop iterates over the bodies *this MPI rank is responsible for*
    for (int i = 0; i < local_N; i++)
    {
        cuerpo1_idx = offset + i; // The actual index of cuerpo1 in the global `cuerpos` array

        // Now, calculate forces between cuerpo1 and *all* other bodies
        for (cuerpo2_idx = 0; cuerpo2_idx < total_N; cuerpo2_idx++)
        {
            if (cuerpo1_idx == cuerpo2_idx)
                continue; // Don't calculate force with self

            if ((all_cuerpos[cuerpo1_idx].px == all_cuerpos[cuerpo2_idx].px) &&
                (all_cuerpos[cuerpo1_idx].py == all_cuerpos[cuerpo2_idx].py) &&
                (all_cuerpos[cuerpo1_idx].pz == all_cuerpos[cuerpo2_idx].pz))
                continue; // Avoid division by zero if bodies are at the same position

            dif_X = all_cuerpos[cuerpo2_idx].px - all_cuerpos[cuerpo1_idx].px;
            dif_Y = all_cuerpos[cuerpo2_idx].py - all_cuerpos[cuerpo1_idx].py;
            dif_Z = all_cuerpos[cuerpo2_idx].pz - all_cuerpos[cuerpo1_idx].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            F = (G * all_cuerpos[cuerpo1_idx].masa * all_cuerpos[cuerpo2_idx].masa) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            // Apply force to cuerpo1_idx. This rank is only calculating the force
            // exerted *on* its assigned bodies (cuerpo1_idx).
            fuerza_totalX[cuerpo1_idx] += dif_X;
            fuerza_totalY[cuerpo1_idx] += dif_Y;
            fuerza_totalZ[cuerpo1_idx] += dif_Z;
        }
    }
}


// This moverCuerpos now operates on a segment of bodies and updates their positions
// based on the already accumulated forces in fuerza_totalX/Y/Z.
// It assumes that fuerza_totalX/Y/Z contain the *final, summed* forces for each body.
void moverCuerpos(cuerpo_t *bodies_to_move, int num_bodies, double delta_tiempo)
{
    for (int i = 0; i < num_bodies; i++)
    {
        // Calculate acceleration
        double inv_masa = (bodies_to_move[i].masa > 0.0f) ? 1.0f / bodies_to_move[i].masa : 0.0f;
        double ax = fuerza_totalX[i] * inv_masa;
        double ay = fuerza_totalY[i] * inv_masa;
        double az = fuerza_totalZ[i] * inv_masa;

        // Update velocities
        bodies_to_move[i].vx += ax * delta_tiempo;
        bodies_to_move[i].vy += ay * delta_tiempo;
        bodies_to_move[i].vz += az * delta_tiempo;

        // Update positions
        bodies_to_move[i].px += bodies_to_move[i].vx * delta_tiempo;
        bodies_to_move[i].py += bodies_to_move[i].vy * delta_tiempo;
        bodies_to_move[i].pz += bodies_to_move[i].vz * delta_tiempo;

        // Reset forces for the next step (important!)
        fuerza_totalX[i] = 0.0;
        fuerza_totalY[i] = 0.0;
        fuerza_totalZ[i] = 0.0;
    }
}


// --- Utility functions (unchanged) ---

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

    cuerpo->r = 1.0;
    cuerpo->g = 1.0;
    cuerpo->b = 1.0;
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

    cuerpo->r = 1.0;
    cuerpo->g = 0.0;
    cuerpo->b = 0.0;
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

    cuerpo->r = 1.0;
    cuerpo->g = 1.0;
    cuerpo->b = 1.0;
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
    // Only the Coordinator frees the main 'cuerpos' array
    if (idW_MPI == 0) {
        free(cuerpos);
    }
}