#include "../inc/nbody.h"

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
}

void moverCuerpos(int idW)
{
    int cuerpo;
    int slice = N / T;
    int ini = idW * slice;
    int lim = ini + slice;

    for (cuerpo = ini; cuerpo < lim; cuerpo++)
    {
        double fx = 0.0, fy = 0.0, fz = 0.0;

        for (int t = 0; t < T; t++)
        {
            fx += fuerza_totalX[t * N + cuerpo];
            fy += fuerza_totalY[t * N + cuerpo];
            fz += fuerza_totalZ[t * N + cuerpo];
        }

        fx /= cuerpos[cuerpo].masa;
        fy /= cuerpos[cuerpo].masa;
        fz /= cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fx * dt;
        cuerpos[cuerpo].vy += fy * dt;
        cuerpos[cuerpo].vz += fz * dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;
        cuerpos[cuerpo].pz += cuerpos[cuerpo].vz * dt;

        for (int t = 0; t < T; t++)
        {
            fuerza_totalX[t * N + cuerpo] = 0.0;
            fuerza_totalY[t * N + cuerpo] = 0.0;
            fuerza_totalZ[t * N + cuerpo] = 0.0;
        }
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
    toroide_r = 5.0;
    toroide_R = 2 * toroide_r;

    srand(0);

    for (int t = 0; t < T; t++)
    {
        for (int i = 0; i < N; i++)
        {
            fuerza_totalX[t * N + i] = 0.0;
            fuerza_totalY[t * N + i] = 0.0;
            fuerza_totalZ[t * N + i] = 0.0;
        }
    }

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

    calcularFuerzas(idW);
    pthread_barrier_wait(&barrera);
    moverCuerpos(idW);
    pthread_barrier_wait(&barrera);

    pthread_exit(NULL);
    return NULL;
}