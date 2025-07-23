#include "../inc/raylib.h"
#include "../inc/nbody.h"

#include <stdio.h>

int main(void)
{
    int screenWidth = GetMonitorWidth(0);
    int screenHeight = GetMonitorHeight(0);
    SetConfigFlags(FLAG_WINDOW_UNDECORATED);
    InitWindow(screenWidth, screenHeight, "N Body Problem - Raylib");

    Camera3D camera = {0};
    camera.position = (Vector3){10.0f, 10.0f, 10.0f};
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    Vector3 origin = {0.0f, 0.0f, 0.0f};

    N = 512;
    pasos = 4096;
    T = 4;

    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (double *)malloc(sizeof(double) * N * T);
    fuerza_totalY = (double *)malloc(sizeof(double) * N * T);
    fuerza_totalZ = (double *)malloc(sizeof(double) * N * T);
    if (cuerpos == NULL || fuerza_totalX == NULL || fuerza_totalY == NULL || fuerza_totalZ == NULL)
    {
        fprintf(stderr, "ERROR: Memory allocation failed\n");
        return 1;
    }
    inicializarCuerpos(cuerpos, N);

    pthread_t threads[T];
    int thread_ids[T];
    pthread_barrier_init(&barrera, NULL, T);

    DisableCursor();

    SetTargetFPS(30);

    while (!WindowShouldClose())
    {
        UpdateCamera(&camera, CAMERA_FREE);

        if (IsKeyPressed('Z'))
            camera.target = (Vector3){0.0f, 0.0f, 0.0f};

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        DrawLine3D(origin, (Vector3){10.0f, 0.0f, 0.0f}, RED);
        DrawLine3D(origin, (Vector3){0.0f, 10.0f, 0.0f}, GREEN);
        DrawLine3D(origin, (Vector3){0.0f, 0.0f, 10.0f}, BLUE);
        DrawGrid(10, 1.0f);

        Vector3 bodypos;
        for (int i = 0; i < N; i++)
        {
            bodypos = (Vector3){cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz};

            switch (cuerpos[i].cuerpo)
            {
            case ESTRELLA:
                DrawSphere(bodypos, 0.1f, RED);
                break;
            case POLVO:
                DrawSphere(bodypos, 0.1f, GREEN);
                break;
            case H2:
                DrawSphere(bodypos, 0.1f, BLUE);
                break;
            default:
                break;
            }
        }

        EndMode3D();

        DrawRectangle(10, 10, 600, 140, Fade(LIME, 0.5f));
        DrawRectangleLines(10, 10, 600, 140, LIME);

        DrawText("Controles por defecto de la cámara libre:", 20, 20, 20, WHITE);
        DrawText("- Rueda del ratón para acercar/alejar", 40, 40, 20, WHITE);
        DrawText("- Presionar la rueda para desplazar", 40, 60, 20, WHITE);
        DrawText("- Z para centrar en (0, 0, 0)", 40, 80, 20, WHITE);
        DrawText("- WASD para mover la cámara", 40, 100, 20, WHITE);
        DrawText("- Mover el ratón para rotar la cámara", 40, 120, 20, WHITE);

        DrawFPS(10, 160);

        dt = GetFrameTime();
        for (int i = 0; i < T; i++)
        {
            thread_ids[i] = i;
            pthread_create(&threads[i], NULL, pfunction, &thread_ids[i]);
        }

        for (int i = 0; i < T; i++)
        {
            pthread_join(threads[i], NULL);
        }

        EndDrawing();
    }

    CloseWindow();

    finalizar();

    return 0;
}