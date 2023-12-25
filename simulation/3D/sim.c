
#include "stdlib.h"
#include "stdio.h"
#include "raylib.h"
#include "raymath.h"

#define SCREEN_WIDTH 512
#define SCREEN_HEIGHT 512

#define c (343.0f*1000)
#define f 40000.0f
#define L (c/f)
#define K 2*PI/L

#define P0 1.0f
#define A 12.0f
#define a 5.0f

#define TRANSDUCER_N 16
#define TRANSDUCER_INTERVAL_XZ (2*a)

#define ARRAY_HEIGHT 250
#define ARRAY_WIDTH (2*TRANSDUCER_N*a)


#define SIM_N_XZ 50
#define SIM_INTERVAL_XZ (ARRAY_WIDTH / SIM_N_XZ)

#define SIM_N_Y 50
#define SIM_INTERVAL_Y (ARRAY_HEIGHT / SIM_N_Y)

Vector3*** createTransducerPts() {
    Vector3*** transducer_pts = (Vector3***)MemAlloc(2*sizeof(Vector3**)); // array for top and bottom

    transducer_pts[0] = (Vector3**)MemAlloc(TRANSDUCER_N*sizeof(Vector3*)); // bottom array
    transducer_pts[1] = (Vector3**)MemAlloc(TRANSDUCER_N*sizeof(Vector3*)); // top array

    for (int i = 0; i < TRANSDUCER_N; i++) {
        transducer_pts[0][i] = (Vector3*)MemAlloc(TRANSDUCER_N*sizeof(Vector3));
        transducer_pts[1][i] = (Vector3*)MemAlloc(TRANSDUCER_N*sizeof(Vector3));
    }

    float x = TRANSDUCER_INTERVAL_XZ/2;
    for (int i = 0; i < TRANSDUCER_N; i++) {
        float z = TRANSDUCER_INTERVAL_XZ/2;
        for (int j = 0; j < TRANSDUCER_N; j++) {
            transducer_pts[0][i][j] = (Vector3){ x, 0.0f, z };
            transducer_pts[1][i][j] = (Vector3){ x, ARRAY_HEIGHT, z };
            z += TRANSDUCER_INTERVAL_XZ;
        }
        x += TRANSDUCER_INTERVAL_XZ;
    }

    return transducer_pts;
}

void freeTransducerPts(Vector3*** transducer_pts) {
    for (int i = 0; i < TRANSDUCER_N; i++) {
        MemFree(transducer_pts[0][i]);
        MemFree(transducer_pts[1][i]);
    }
    MemFree(transducer_pts[0]);
    MemFree(transducer_pts[1]);
    MemFree(transducer_pts);
}

Vector3*** createSimPts() {
    Vector3*** sim_pts = (Vector3***)MemAlloc(SIM_N_Y*sizeof(Vector3**)); // array for each vertical slice
    for (int i = 0; i < SIM_N_Y; i++) {
        sim_pts[i] = (Vector3**)MemAlloc(SIM_N_XZ*sizeof(Vector3*));
        for (int j = 0; j < SIM_N_XZ; j++) {
            sim_pts[i][j] = (Vector3*)MemAlloc(SIM_N_XZ*sizeof(Vector3));
        }
    }

    float y = SIM_INTERVAL_Y/2;
    for (int i = 0; i < SIM_N_Y; i++) {
        float x = SIM_INTERVAL_XZ/2;
        for (int j = 0; j < SIM_N_XZ; j++) {
            float z = SIM_INTERVAL_XZ/2;
            for (int k = 0; k < SIM_N_XZ; k++) {
                Vector3 pt = { x, y, z };
                sim_pts[i][j][k] = pt;
                z += SIM_INTERVAL_XZ;
            }
            x += SIM_INTERVAL_XZ;
        }
        y += SIM_INTERVAL_Y;
    }

    return sim_pts;
}

float*** createPowerPts() {
    float*** power_pts = (float***)MemAlloc(SIM_N_Y*sizeof(float**)); // array for each vertical slice
    for (int i = 0; i < SIM_N_Y; i++) {
        power_pts[i] = (float**)MemAlloc(SIM_N_XZ*sizeof(float*));
        for (int j = 0; j < SIM_N_XZ; j++) {
            power_pts[i][j] = (float*)MemAlloc(SIM_N_XZ*sizeof(float));
        }
    }

    for (int i = 0; i < SIM_N_Y; i++) {
        for (int j = 0; j < SIM_N_XZ; j++) {
            for (int k = 0; k < SIM_N_XZ; k++) {
                power_pts[i][j][k] = 0.0f;
            }
        }
    }

    return power_pts;
}

float*** createPhiPts() {
    float*** phi_pts = (float***)MemAlloc(2*sizeof(float**)); // array for top and bottom

    phi_pts[0] = (float**)MemAlloc(TRANSDUCER_N*sizeof(float*)); // bottom array
    phi_pts[1] = (float**)MemAlloc(TRANSDUCER_N*sizeof(float*)); // top array

    for (int i = 0; i < TRANSDUCER_N; i++) {
        phi_pts[0][i] = (float*)MemAlloc(TRANSDUCER_N*sizeof(float));
        phi_pts[1][i] = (float*)MemAlloc(TRANSDUCER_N*sizeof(float));
    }

    for (int i = 0; i < TRANSDUCER_N; i++) {
        for (int j = 0; j < TRANSDUCER_N; j++) {
            phi_pts[0][i][j] = 0.0f;
            phi_pts[1][i][j] = 0.0f;
        }
    }

    return phi_pts;
}

void freePhiPts(float*** phi_pts) {
    for (int i = 0; i < TRANSDUCER_N; i++) {
        MemFree(phi_pts[0][i]);
        MemFree(phi_pts[1][i]);
    }
    MemFree(phi_pts[0]);
    MemFree(phi_pts[1]);
    MemFree(phi_pts); 
}

void freePowerPts(float*** power_pts) {

    for (int i = 0; i < SIM_N_Y; i++) {
        for (int j = 0; j < SIM_N_XZ; j++) {
            MemFree(power_pts[i][j]);
        }
        MemFree(power_pts[i]);
    }

    MemFree(power_pts);    
}

void freeSimPts(Vector3*** sim_pts) {

    for (int i = 0; i < SIM_N_Y; i++) {
        for (int j = 0; j < SIM_N_XZ; j++) {
            MemFree(sim_pts[i][j]);
        }
        MemFree(sim_pts[i]);
    }

    MemFree(sim_pts);    
}

float sinc(float x) {
    if (x == 0) return 1.0f;
    return sinf(PI*x) / (PI*x);
}

float Df(float theta) {
    return sinc(K*a*sinf(theta));
}

float calcTheta(Vector3 v, Vector3 u) {
    Vector3 v_norm = Vector3Normalize(v);
    Vector3 u_norm = Vector3Normalize(u);
    float dot = Vector3DotProduct(v_norm, u_norm);

    return acosf(dot);
}

int P(Vector3 pt, Vector3 t_pt, bool isTop, float t_phi, Vector2* complex) {

    Vector3 t_norm = (isTop) ? (Vector3){ 0.0f, 1.0f, 0.0f } : (Vector3){ 0.0f, -1.0f, 0.0f };

    float d = Vector3Distance(pt, t_pt);
    if (d == 0) {
        return 0;
    }
    
    Vector3 v1 = Vector3Add(pt, t_pt);
    Vector3 v2 = Vector3Add(t_norm, t_pt);
    float theta = calcTheta(v1, v2);

    float real = P0*A*Df(theta)*cosf(t_phi + K*d);
    float im = P0*A*Df(theta)*sinf(t_phi + K*d);

    complex->x += real;
    complex->y += im;

    return 0;
}

float calculatePhase(Vector3 focus, Vector3 t_pt, bool isTop) {

    float d = Vector3Distance(focus, t_pt);
    float phi = -d*K;
    
    if (isTop) phi += PI;

    // phi *= (32 / PI);
    // phi = fmodf(phi, 64);
    // printf("phi: %f", phi);

    phi = fmodf(phi, 2*PI);
    return phi;
}

int calculatePower(float*** power,
                   Vector3*** sim_pts,
                   Vector3*** transducer_pts,
                   float*** phi_pts) {
    float max = -1.0f;
    float min = INFINITY;

    Vector2 *complex = calloc(1, sizeof(Vector2));

    for (int i = 0; i < SIM_N_Y; i++) {
        for (int j = 0; j < SIM_N_XZ; j++) {
            for (int k = 0; k < SIM_N_XZ; k++) {

                complex->x = 0.0f;
                complex->y = 0.0f;

                Vector3 pt = sim_pts[i][j][k];

                for (int m = 0; m < TRANSDUCER_N; m++) {
                    for (int n = 0; n < TRANSDUCER_N; n++) {
                        Vector3 t_pt = transducer_pts[0][m][n];
                        P(pt, t_pt, false, phi_pts[0][m][n], complex);

                        t_pt = transducer_pts[1][m][n];
                        P(pt, t_pt, true, phi_pts[1][m][n], complex);
                    }
                }

                float p = Vector2Length(*complex);

                if (p > max) {
                    max = p;
                }
                if (p < min) {
                    min = p;
                }
                power[i][j][k] = p;
            }
        }
    }

    for (int i = 0; i < SIM_N_Y; i++) {
        for (int j = 0; j < SIM_N_XZ; j++) {
            for (int k = 0; k < SIM_N_XZ; k++) {
                float p = power[i][j][k];
                p -= min;
                p /= max;
                p *= 255.0f;
                p = roundf(p);

                power[i][j][k] = p;
            }
        }
    }

    free(complex);

    return 0;
}

int main() {

    Vector3*** transducer_pts = createTransducerPts();
    Vector3*** sim_pts = createSimPts();
    float*** power_pts = createPowerPts();
    float*** phi_pts = createPhiPts();


    Vector3 focus = {ARRAY_WIDTH/2, ARRAY_HEIGHT/2, ARRAY_WIDTH/2};

    for (int i = 0; i < TRANSDUCER_N; i++) {
        for (int j = 0; j < TRANSDUCER_N; j++) {
            Vector3 t_pt = transducer_pts[0][i][j];
            float phi = calculatePhase(focus, t_pt, false);
            phi_pts[0][i][j] = phi;

            t_pt = transducer_pts[1][i][j];
            phi = calculatePhase(focus, t_pt, true);
            phi_pts[1][i][j] = phi;
        }
    }

    calculatePower(power_pts, sim_pts, transducer_pts, phi_pts);

    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "3D Acoustic Simulation");

    Camera camera = { 0 };
    camera.position = (Vector3){ 3*ARRAY_WIDTH, ARRAY_HEIGHT/2, 3*ARRAY_WIDTH };
    camera.target = (Vector3){ ARRAY_WIDTH/2, ARRAY_HEIGHT/2, ARRAY_WIDTH/2 };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    SetTargetFPS(30);

    bool exit = false;
    while (!exit && !WindowShouldClose()) {

        // UpdateCamera(&camera, CAMERA_ORBITAL);

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode3D(camera);

        /*
            Y
            |___ X
           /
          Z  
        */

        for (int i = 0; i < TRANSDUCER_N; i++) {
            for (int j = 0; j < TRANSDUCER_N; j++) {
                Vector3 pt = transducer_pts[0][i][j];

                DrawCylinder(pt, a, a, -a/2, 25, GREEN);
                DrawCylinderWires(pt, a, a, -a/2, 25, DARKGREEN);

                pt = transducer_pts[1][i][j];
                DrawCylinder(pt, a, a, a/2, 25, GREEN);
                DrawCylinderWires(pt, a, a, a/2, 25, DARKGREEN);
            }
        }

        for (int i = 0; i < SIM_N_Y; i++) {
            for (int j = 0; j < SIM_N_XZ; j++) {
                for (int k = 0; k < SIM_N_XZ; k++) {
                    Vector3 pt = sim_pts[i][j][k];

                    // DrawCube(pt, SIM_INTERVAL_XZ, SIM_INTERVAL_Y, SIM_INTERVAL_XZ, PURPLE);
                    // DrawCubeWires(pt, SIM_INTERVAL_XZ, SIM_INTERVAL_Y, SIM_INTERVAL_XZ, DARKPURPLE);

                    float p = power_pts[i][j][k];
                    if (p > 50) {
                        Color color = {p, p, p, 0x80}; // ColorFromHSV(p, 100.0, 100.0);
                        DrawSphere(pt, 1.5f, color);
                    }
                }
            }
        }
        DrawSphere(focus, 1.0f, RED);

        EndMode3D();
        EndDrawing();
        
        // TakeScreenshot(const char *fileName); 

        // exit = true;

    }

    CloseWindow();

    freeTransducerPts(transducer_pts);
    freeSimPts(sim_pts);
    freePowerPts(power_pts);
    freePhiPts(phi_pts);

    return 0;
}