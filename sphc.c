#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"

// Current time
int steps;
double time;
double dt;

// Particle properties
// Number of particles
int n;

// Original position
double * ox;
double * oy;
double * oz;

// Position
double * px;
double * py;
double * pz;

// Velocity
double * vx;
double * vy;
double * vz;

// Force
double * fx;
double * fy;
double * fz;

// Neighbors
int gridsize;
int * gridc;
int ** grid;
int * nc;
int ** nbs;

// Candidate position
double * cpx;
double * cpy;
double * cpz;

// Candidate velocity
double * cvx;
double * cvy;
double * cvz;

// Lambda
double * lm;

// Pressure
double * dpx;
double * dpy;
double * dpz;

// Vorticity force
double * vox;
double * voy;
double * voz;

// Write output
void write_step(int step) {
    char fname[200];
    sprintf(fname, "output/particles_%05d", step);
    FILE * fout = fopen(fname, "w");
    fprintf(fout, "%d\n", n);
    for (int i = 0; i < n; ++i) {
        fprintf(fout, "%lf %lf %lf %lf %lf %lf\n",
            px[i], py[i], pz[i], vx[i], vy[i], vz[i]);
    }
    fclose(fout);
}

void write_candidates(int step) {
    char fname[200];
    sprintf(fname, "output/candidates_%05d", step);
    FILE * fout = fopen(fname, "w");
    fprintf(fout, "%d\n", n);
    for (int i = 0; i < n; ++i) {
        fprintf(fout, "%lf %lf %lf %lf %lf %lf\n",
            cpx[i], cpy[i], cpz[i], cvx[i], cvy[i], cvz[i]);
    }
    fclose(fout);
}

// Runtime computed constant
double PRESSURE_RADIUS_FACTOR = 0.0;

void computePressureRadiusFactor() {
    const double c = 315.0 / (64.0 * PI);
    double r = ARTIFICIAL_PRESSURE_RADIUS;
    PRESSURE_RADIUS_FACTOR = 1.0 / (c * pow(
        KERNEL_SIZE * KERNEL_SIZE - r * r, 3
    ) / pow(KERNEL_SIZE, 9));
}

// 6th degree polynomial kernel
double poly6kernel(int i, int j) {
    const double c = 315.0 / (64.0 * PI);
    double dx = cpx[i] - cpx[j];
    double dy = cpy[i] - cpy[j];
    double dz = cpz[i] - cpz[j];
    double r = sqrt(dx * dx + dy * dy + dz * dz);

    return (r > KERNEL_SIZE) ? 0 : c * pow(
        KERNEL_SIZE * KERNEL_SIZE - r * r, 3
    ) / pow(KERNEL_SIZE, 9);
}

// Spiky gradient kernel
void spikykernel(int i, int j, double gv[]) {
    // IS THIS THE WRONG C Constant? Should be 3*15 / PI? (+corrected)
    // This was for spiky kernel in Java code, but this is actually spikyGrad?
    const double c = 3.0 * 15.0 / PI;
    double dx = cpx[i] - cpx[j];
    double dy = cpy[i] - cpy[j];
    double dz = cpz[i] - cpz[j];
    double r = sqrt(dx * dx + dy * dy + dz * dz);
    if (r < 1e-4) {
        r = 1e-4;
    }

    if (r > KERNEL_SIZE) {
        gv[0] = 0.0;
        gv[1] = 0.0;
        gv[2] = 0.0;
    }
    else {
        double f = c / pow(KERNEL_SIZE, 6) * pow(KERNEL_SIZE - r, 2) / r;
        gv[0] = f * dx;
        gv[1] = f * dy;
        gv[2] = f * dz;
    }
}

// Viscosity kernel
double viscositykernel(int i, int j) {
    // IS THIS THE WRONG C Constant? Should be 15.0 / (2*PI)?
    // (- java code was sort of messed up, will play around with this later)
    const double c = 315.0 / (64.0 * PI);
    double dx = cvx[i] - cvx[j];
    double dy = cvy[i] - cvy[j];
    double dz = cvz[i] - cvz[j];
    double r = sqrt(dx * dx + dy * dy + dz * dz);
  
    // THIS ISNT RIGHT. THIS IS FROM POLY6
    return (r > KERNEL_SIZE) ? 0 : c * pow(
        KERNEL_SIZE * KERNEL_SIZE - r * r, 3
    ) / pow(KERNEL_SIZE, 9);
}

// Advance particles by a single step
void step() {
    // Clear forces
    memset(fx, 0, n * sizeof(double));
    memset(fy, 0, n * sizeof(double));
    memset(fz, 0, n * sizeof(double));

    // Apply forces
    for (int i = 0; i < n; ++i) {
        // Gravity
        fy[i] += GRAVITY;

        // Vorticity
        fx[i] += vox[i];
        fy[i] += voy[i];
        fz[i] += voz[i];
    }

    // Compute candidate velocities and positions
    for (int i = 0; i < n; ++i) {
        cvx[i] = vx[i] + fx[i] * dt;
        cvy[i] = vy[i] + fy[i] * dt;
        cvz[i] = vz[i] + fz[i] * dt;

        // TODO Velocity dampening (if used) goes here

        cpx[i] = px[i] + cvx[i] * dt;
        cpy[i] = py[i] + cvy[i] * dt;
        cpz[i] = pz[i] + cvz[i] * dt;
    }

    // Find neighbors
    // Place into grid
    memset(gridc, 0, n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        int gx = (px[i] / KERNEL_SIZE);
        int gy = (py[i] / KERNEL_SIZE);
        int gz = (pz[i] / KERNEL_SIZE);
        if (gx < 0) gx = 0; else if (gx >= gridsize) gx = gridsize - 1;
        if (gy < 0) gy = 0; else if (gy >= gridsize) gy = gridsize - 1;
        if (gz < 0) gz = 0; else if (gz >= gridsize) gz = gridsize - 1;

        int offset = gx * gridsize * gridsize + gy * gridsize + gz;
        grid[offset][gridc[offset]] = i;
        gridc[offset]++;
    }

    // Find neighbors using the grid
    memset(nc, 0, n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        int gx = (px[i] / KERNEL_SIZE);
        int gy = (py[i] / KERNEL_SIZE);
        int gz = (pz[i] / KERNEL_SIZE);

        // Look at surrounding grid spaces
        for (int gox = gx - 1; gox <= gx + 1; ++gox) {
        for (int goy = gy - 1; goy <= gy + 1; ++goy) {
        for (int goz = gz - 1; goz <= gz + 1; ++goz) {
            if (gox < 0 || gox >= gridsize) continue;
            if (goy < 0 || goy >= gridsize) continue;
            if (goz < 0 || goz >= gridsize) continue;
            int offset = gox * gridsize * gridsize + goy * gridsize + goz;
            // Iterate through grid space
            for (int j = 0, l = gridc[offset]; j < l; ++j) {
                // Check distance
                int nb = grid[offset][j];
                double dx = cpx[i] - cpx[nb];
                double dy = cpy[i] - cpy[nb];
                double dz = cpz[i] - cpz[nb];
                double mag = sqrt(dx * dx + dy * dy + dz * dz);
                if (nb != i && mag <= KERNEL_SIZE) {
                    // Add as neighbor
                    nbs[i][nc[i]] = nb;
                    nc[i]++;
                }
            }
        }
        }
        }
    }

    for (int iteration = 0; iteration < ITERATIONS_PER_STEP; iteration++) {
        // Compute lambda
        for (int i = 0; i < n; ++i) {
            double rho = 0.0;
            double numer = 0.0;
            double denom = 0.0;

            double sx = 0.0;
            double sy = 0.0;
            double sz = 0.0;
            double gv[3];

            for (int j = 0, l = nc[i]; j < l; ++j) {
                int nb = nbs[i][j];
                rho += PARTICLE_MASS * poly6kernel(i, nb);

                spikykernel(i, nb, gv);
                sx += gv[0];
                sy += gv[1];
                sz += gv[2];

                gv[0] /= -PARTICLE_DENSITY;
                gv[1] /= -PARTICLE_DENSITY;
                gv[2] /= -PARTICLE_DENSITY;

                denom += gv[0] * gv[0] + gv[1] * gv[1] + gv[2] * gv[2];
            }
            numer = rho / PARTICLE_DENSITY - 1.0;
            sx /= PARTICLE_DENSITY;
            sy /= PARTICLE_DENSITY;
            sz /= PARTICLE_DENSITY;

            denom += sx * sx + sy * sy + sz * sz;

            lm[i] = -numer / (denom + FORCE_EPSILON);
        }

        // Jiggle particles
        for (int i = 0; i < n; ++i) {
            // Compute pressure
            double sx = 0.0;
            double sy = 0.0;
            double sz = 0.0;
            double gv[3];
            double c;
            for (int j = 0, l = nc[i]; j < l; ++j) {
                int nb = nbs[i][j];
                double scorr = 0.0;
                // Surface tension
                // TODO Enable / disable surface tension
                scorr = -ARTIFICIAL_PRESSURE_STRENGTH * pow(
                    poly6kernel(i, nb) * PRESSURE_RADIUS_FACTOR,
                    ARTIFICIAL_PRESSURE_POWER
                );

                spikykernel(i, nb, gv);
                c = lm[i] + lm[nb] + scorr;
                sx += gv[0] * c;
                sy += gv[1] * c;
                sz += gv[2] * c;
            }

            dpx[i] = sx / PARTICLE_DENSITY;
            dpy[i] = sy / PARTICLE_DENSITY;
            dpz[i] = sz / PARTICLE_DENSITY;

            // Update candidate position
            cpx[i] += dpx[i];
            cpy[i] += dpy[i];
            cpz[i] += dpz[i];

            cvx[i] = (cpx[i] - px[i]) / dt;
            cvy[i] = (cpy[i] - py[i]) / dt;
            cvz[i] = (cpz[i] - pz[i]) / dt;

            // Resolve collisions
            if (cpy[i] < 0.0) { // Floor
                cpy[i] = 0.0;
                cvy[i] = -cvy[i];
            }

            if (cpx[i] < 0.0) { // Wall
                cpx[i] = 0.0;
                cvx[i] = -cvx[i];
            }
            if (cpx[i] > BOX_SIZE) { // Wall
                cpx[i] = BOX_SIZE;
                cvx[i] = -cvx[i];
            }
            if (cpz[i] < 0.0) { // Wall
                cpz[i] = 0.0;
                cvz[i] = -cvz[i];
            }
            if (cpz[i] > BOX_SIZE) { // Wall
                cpz[i] = BOX_SIZE;
                cvz[i] = -cvz[i];
            }
        }
    }

    // Compute angular velocity
    for (int i = 0; i < n; ++i) {
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        double gv[3];
        double velx, vely, velz;
        for (int j = 0, l = nc[i]; j < l; ++j) {
            int nb = nbs[i][j];
          // I THINK VX, VY, VZ SHOULDNT BE VELOCITY BUT A DIFFERENT VARIABLE FOR VORTICITY
          // (+corrected)
            velx = cvx[nb] - cvx[i];
            vely = cvy[nb] - cvy[i];
            velz = cvz[nb] - cvz[i];

            spikykernel(i, nb, gv);
          // SAME QUESTION HERE FOR VX VY VZ
            sx += vely * gv[2] - velz * gv[1];
            sy += velz * gv[0] - velx * gv[2];
            sz += velx * gv[1] - vely * gv[0];
        }

        // Reusing dp array
        dpx[i] = sx;
        dpy[i] = sy;
        dpz[i] = sz;
    }

    // Compute vorticity force
    for (int i = 0; i < n; ++i) {
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        double gv[3];
        double mag;
        for (int j = 0, l = nc[i]; j < l; ++j) {
            int nb = nbs[i][j];
            spikykernel(i, nb, gv);
            // USED L TWICE IN FOR LOOP AND IN BELOW (+corrected)
            mag = sqrt(dpx[nb] * dpx[nb] + dpy[nb] * dpy[nb] + dpz[nb] * dpz[nb]);
            gv[0] *= mag / PARTICLE_DENSITY;
            gv[1] *= mag / PARTICLE_DENSITY;
            gv[2] *= mag / PARTICLE_DENSITY;

            sx += gv[0];
            sy += gv[1];
            sz += gv[2];
        }
        mag = sqrt(sx * sx + sy * sy + sz * sz);
        if (mag < 1e-5) {
            mag = 1e-5;
        }
        sx /= mag;
        sy /= mag;
        sz /= mag;

        vox[i] = (sy * dpz[i] - sz * dpy[i]) * VORTICITY_COEFFICIENT;
        voy[i] = (sz * dpx[i] - sx * dpz[i]) * VORTICITY_COEFFICIENT;
        // I THINK BELOW SHOULD BE SX*DPY - SY*DPX (+corrected)
        voz[i] = (sx * dpy[i] - sy * dpx[i]) * VORTICITY_COEFFICIENT;
    }

    // Viscosity
    for (int i = 0; i < n; ++i) {
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        double dx, dy, dz;
        for (int j = 0, l = nc[i]; j < l; ++j) {
            int nb = nbs[i][j];
            dx = cvx[nb] - cvx[i];
            dy = cvy[nb] - cvy[i];
            dz = cvz[nb] - cvz[i];
            // I THINK THE POLY6KERNEL IS USED IN THE JAVA VERSION?
            double c = viscositykernel(i, nb);
            sx += c * dx;
            sy += c * dy;
            sz += c * dz;
        }
        cvx[i] += ARTIFICIAL_VISCOSITY * sx;
        cvy[i] += ARTIFICIAL_VISCOSITY * sy;
        cvz[i] += ARTIFICIAL_VISCOSITY * sz;
    }

    // Update particle
    for (int i = 0; i < n; ++i) {
        vx[i] = cvx[i];
        vy[i] = cvy[i];
        vz[i] = cvz[i];

        px[i] = cpx[i];
        py[i] = cpy[i];
        pz[i] = cpz[i];
    }
}

void run() {
    time = 0.0;
    write_step(0);
    for (int i = 1; i <= steps; ++i) {
        step();
        time += dt;
        write_step(i);
    }
}

void init(const char * initfile) {
    FILE * fin = fopen(initfile, "r");
    // Read number of particles
    fscanf(fin, "%d", &n);

    // Allocate arrays
    ox = (double *) malloc(n * sizeof(double));
    oy = (double *) malloc(n * sizeof(double));
    oz = (double *) malloc(n * sizeof(double));

    px = (double *) malloc(n * sizeof(double));
    py = (double *) malloc(n * sizeof(double));
    pz = (double *) malloc(n * sizeof(double));
    
    vx = (double *) malloc(n * sizeof(double));
    vy = (double *) malloc(n * sizeof(double));
    vz = (double *) malloc(n * sizeof(double));

    fx = (double *) malloc(n * sizeof(double));
    fy = (double *) malloc(n * sizeof(double));
    fz = (double *) malloc(n * sizeof(double));

    gridsize = (int) (BOX_SIZE / KERNEL_SIZE + 0.5);
    gridc = (int *) malloc(gridsize * gridsize * gridsize * n * sizeof(int));
    grid = (int **) malloc(gridsize * gridsize * gridsize * sizeof(int *));
    for (int i = 0; i < gridsize * gridsize * gridsize; ++i) {
        grid[i] = (int *) malloc(n * sizeof(int));
    }
    nc = (int *) malloc(n * sizeof(int));
    nbs = (int **) malloc(n * sizeof(int *));
    for (int i = 0; i < n; ++i) {
        nbs[i] = (int *) malloc(n * sizeof(int));
    }

    cpx = (double *) malloc(n * sizeof(double));
    cpy = (double *) malloc(n * sizeof(double));
    cpz = (double *) malloc(n * sizeof(double));

    cvx = (double *) malloc(n * sizeof(double));
    cvy = (double *) malloc(n * sizeof(double));
    cvz = (double *) malloc(n * sizeof(double));

    lm = (double *) malloc(n * sizeof(double));

    dpx = (double *) malloc(n * sizeof(double));
    dpy = (double *) malloc(n * sizeof(double));
    dpz = (double *) malloc(n * sizeof(double));

    vox = (double *) malloc(n * sizeof(double));
    voy = (double *) malloc(n * sizeof(double));
    voz = (double *) malloc(n * sizeof(double));
    memset(vox, 0, n * sizeof(double));
    memset(voy, 0, n * sizeof(double));
    memset(voz, 0, n * sizeof(double));

    // Read initial positions
    for (int i = 0; i < n; ++i) {
        fscanf(fin, "%lf %lf %lf %lf %lf %lf",
            px + i, py + i, pz + i, vx + i, vy + i, vz + i);
        ox[i] = px[i];
        oy[i] = py[i];
        oz[i] = pz[i];
    }

    // Compute pressure radius factor
    computePressureRadiusFactor();

    fclose(fin);
}

int main(int argc, char ** argv) {
    if (argc < 4) {
        printf("./sphc <init file> <steps> <dt>\n");
        return 1;
    }

    init(argv[1]);
    steps = atoi(argv[2]);
    dt = atof(argv[3]);
    run();

	return 0;
}

