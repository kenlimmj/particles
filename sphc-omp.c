#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "constants.h"

// Current time
int steps;

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
int * gx;
int * gy;
int * gz;
omp_lock_t * gridlock;
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

// Timing
double time_force = 0.0;
double time_candidate = 0.0;
double time_neighbor = 0.0;
double time_solve = 0.0;
double time_vorticity = 0.0;
double time_viscosity = 0.0;
double time_update = 0.0;
double total_time_force = 0.0;
double total_time_candidate = 0.0;
double total_time_neighbor = 0.0;
double total_time_solve = 0.0;
double total_time_vorticity = 0.0;
double total_time_viscosity = 0.0;
double total_time_update = 0.0;

// Write output
void write_step(int step) {
    char fname[200];
    sprintf(fname, "output/data/particles_%05d", step);
    FILE * fout = fopen(fname, "w");
    fprintf(fout, "%d\n", n);
    fprintf(fout, "%lf\n", BOX_SIZE);
    for (int i = 0; i < n; ++i) {
        fprintf(fout, "%lf %lf %lf %lf %lf %lf\n",
            px[i], py[i], pz[i], vx[i], vy[i], vz[i]);
    }
    fclose(fout);
}

void write_candidates(int step) {
    char fname[200];
    sprintf(fname, "output/data/candidates_%05d", step);
    FILE * fout = fopen(fname, "w");
    fprintf(fout, "%d\n", n);
    fprintf(fout, "%lf\n", BOX_SIZE);
    for (int i = 0; i < n; ++i) {
        fprintf(fout, "%lf %lf %lf %lf %lf %lf\n",
            cpx[i], cpy[i], cpz[i], cvx[i], cvy[i], cvz[i]);
    }
    fclose(fout);
}

// Runtime computed constant
double PRESSURE_RADIUS_FACTOR = 0.0;

void computeConstants() {
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
void viscositykernel(int i, int j, double gv[]) {
    const double c = 315.0 / (64.0 * PI);
    double dx = cvx[i] - cvx[j];
    double dy = cvy[i] - cvy[j];
    double dz = cvz[i] - cvz[j];
    double r = sqrt(dx * dx + dy * dy + dz * dz);
    
    double k = (r > KERNEL_SIZE) ? 0 : c * pow(
        KERNEL_SIZE * KERNEL_SIZE - r * r, 3
    ) / pow(KERNEL_SIZE, 9);
    
    gv[0] = k * dx;
    gv[1] = k * dy;
    gv[2] = k * dz;
}

// Advance particles by a single step
void step(double dt) {
    double start;
    #pragma omp parallel shared (start)
    {
    
    #pragma omp single nowait
    start = omp_get_wtime();

    // Apply forces
    #pragma omp for
    for (int i = 0; i < n; ++i) {
        // Gravity
        fx[i] = 0.0;
        fy[i] = GRAVITY;
        fz[i] = 0.0;

        // Vorticity
        fx[i] = vox[i];
        fy[i] += voy[i];
        fz[i] = voz[i];
    }
    
    #pragma omp single nowait
    {
    time_force = omp_get_wtime() - start;
    start = omp_get_wtime();
    }
    
    // Compute candidate velocities and positions
    #pragma omp for
    for (int i = 0; i < n; ++i) {
        cvx[i] = vx[i] + fx[i] * dt;
        cvy[i] = vy[i] + fy[i] * dt;
        cvz[i] = vz[i] + fz[i] * dt;

        // TODO Velocity dampening (if used) goes here

        cpx[i] = px[i] + cvx[i] * dt;
        cpy[i] = py[i] + cvy[i] * dt;
        cpz[i] = pz[i] + cvz[i] * dt;
    }
    
    #pragma omp single nowait
    {
    time_candidate = omp_get_wtime() - start;
    start = omp_get_wtime();
    }
    
    // Find neighbors
    // Place into grid
    #pragma omp single
    memset(gridc, 0, gridsize * gridsize * gridsize * sizeof(int));
    #pragma omp for
    for (int i = 0; i < n; ++i) {
        gx[i] = (px[i] / KERNEL_SIZE);
        gy[i] = (py[i] / KERNEL_SIZE);
        gz[i] = (pz[i] / KERNEL_SIZE);
        if (gx[i] < 0) gx[i] = 0; else if (gx[i] >= gridsize) gx[i] = gridsize - 1;
        if (gy[i] < 0) gy[i] = 0; else if (gy[i] >= gridsize) gy[i] = gridsize - 1;
        if (gz[i] < 0) gz[i] = 0; else if (gz[i] >= gridsize) gz[i] = gridsize - 1;
        
        int offset = gx[i] * gridsize * gridsize + gy[i] * gridsize + gz[i];
        omp_set_lock(&gridlock[offset]);
        grid[offset][gridc[offset]] = i;
        gridc[offset]++;
        omp_unset_lock(&gridlock[offset]);
    }

    // Find neighbors using the grid
    #pragma omp single
    memset(nc, 0, n * sizeof(int));
    #pragma omp for
    for (int i = 0; i < n; ++i) {
        // Look at surrounding grid spaces
        for (int gox = gx[i] - 1; gox <= gx[i] + 1; ++gox) {
        for (int goy = gy[i] - 1; goy <= gy[i] + 1; ++goy) {
        for (int goz = gz[i] - 1; goz <= gz[i] + 1; ++goz) {
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
    
    #pragma omp single nowait
    {
    time_neighbor = omp_get_wtime() - start;
    start = omp_get_wtime();
    }
    
    for (int iteration = 0; iteration < ITERATIONS_PER_STEP; iteration++) {
        // Compute lambda
        #pragma omp for
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
        #pragma omp for
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
                cvy[i] = -0.75 * cvy[i];
            }

            if (cpx[i] < 0.0) { // Wall
                cpx[i] = 0.0;
                cvx[i] = -0.75 * cvx[i];
            }
            if (cpx[i] > BOX_SIZE) { // Wall
                cpx[i] = BOX_SIZE;
                cvx[i] = -0.75 * cvx[i];
            }
            if (cpz[i] < 0.0) { // Wall
                cpz[i] = 0.0;
                cvz[i] = -0.75 * cvz[i];
            }
            if (cpz[i] > BOX_SIZE) { // Wall
                cpz[i] = BOX_SIZE;
                cvz[i] = -0.75 * cvz[i];
            }
        }
    }
    
    #pragma omp single nowait
    {
    time_solve = omp_get_wtime() - start;
    start = omp_get_wtime();
    }
    
    // Compute angular velocity
    #pragma omp for
    for (int i = 0; i < n; ++i) {
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        double gv[3];
        double velx, vely, velz;
        for (int j = 0, l = nc[i]; j < l; ++j) {
            int nb = nbs[i][j];
            velx = cvx[nb] - cvx[i];
            vely = cvy[nb] - cvy[i];
            velz = cvz[nb] - cvz[i];

            spikykernel(i, nb, gv);
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
    #pragma omp for
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
        voz[i] = (sx * dpy[i] - sy * dpx[i]) * VORTICITY_COEFFICIENT;
    }
    
    #pragma omp single nowait
    {
    time_vorticity = omp_get_wtime() - start;
    start = omp_get_wtime();
    }
    
    // Viscosity
    #pragma omp for
    for (int i = 0; i < n; ++i) {
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        double dx, dy, dz;
        double gv[3];
        for (int j = 0, l = nc[i]; j < l; ++j) {
            int nb = nbs[i][j];
            viscositykernel(i, nb, gv);
            sx += gv[0];
            sy += gv[1];
            sz += gv[2];
        }
        cvx[i] += ARTIFICIAL_VISCOSITY * sx;
        cvy[i] += ARTIFICIAL_VISCOSITY * sy;
        cvz[i] += ARTIFICIAL_VISCOSITY * sz;
    }
    
    #pragma omp single nowait
    {
    time_viscosity = omp_get_wtime() - start;
    start = omp_get_wtime();
    }
    
    // Update particle
    #pragma omp for
    for (int i = 0; i < n; ++i) {
        vx[i] = cvx[i];
        vy[i] = cvy[i];
        vz[i] = cvz[i];

        px[i] = cpx[i];
        py[i] = cpy[i];
        pz[i] = cpz[i];
    }
    
    #pragma omp single nowait
    {
    time_update = omp_get_wtime() - start;
    
    /*
    printf("Forces:     %lf\n", time_force);
    printf("Candidates: %lf\n", time_candidate);
    printf("Neighbors:  %lf\n", time_neighbor);
    printf("Iterations: %lf\n", time_solve);
    printf("Vorticity:  %lf\n", time_vorticity);
    printf("Viscosity:  %lf\n", time_viscosity);
    printf("Update:     %lf\n", time_update);
    printf("\n");
    */
    
    total_time_force += time_force;
    total_time_candidate += time_candidate;
    total_time_neighbor += time_neighbor;
    total_time_solve += time_solve;
    total_time_vorticity += time_vorticity;
    total_time_viscosity += time_viscosity;
    total_time_update += time_update;
    }
    }
}

void run(int steps, double dt) {
    double time = 0.0;
    write_step(0);
    double sdt = dt / (double) STEPS_PER_FRAME;
    for (int i = 1; i <= steps; ++i) {
        for (int j = 0; j < STEPS_PER_FRAME; ++j) {
            step(sdt);
            time += sdt;
        }
        write_step(i);
    }
    
    printf("Total Forces:     %lf\n", total_time_force);
    printf("Total Candidates: %lf\n", total_time_candidate);
    printf("Total Neighbors:  %lf\n", total_time_neighbor);
    printf("Total Iterations: %lf\n", total_time_solve);
    printf("Total Vorticity:  %lf\n", total_time_vorticity);
    printf("Total Viscosity:  %lf\n", total_time_viscosity);
    printf("Total Update:     %lf\n", total_time_update);
    printf("Total Time: %lf\n", total_time_force +
        total_time_candidate + total_time_neighbor + total_time_solve +
        total_time_vorticity + total_time_viscosity + total_time_update);
    printf("\n");
}

void init(const char * initfile) {
    FILE * fin = fopen(initfile, "r");
    // Read number of particles
    fscanf(fin, "%d", &n);
    fscanf(fin, "%lf", &BOX_SIZE);

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
    gx = (int *) malloc(n * sizeof(int));
    gy = (int *) malloc(n * sizeof(int));
    gz = (int *) malloc(n * sizeof(int));
    gridlock = (omp_lock_t *) malloc(gridsize * gridsize * gridsize * sizeof(omp_lock_t));
    for (int i = 0; i < gridsize * gridsize * gridsize; ++i) {
        omp_init_lock(&gridlock[i]);
    }
    gridc = (int *) malloc(gridsize * gridsize * gridsize * sizeof(int));
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

    computeConstants();

    fclose(fin);
}

int main(int argc, char ** argv) {
    if (argc < 4) {
        printf("./sphc <init file> <steps> <dt>\n");
        return 1;
    }

    init(argv[1]);
    double steps = atoi(argv[2]);
    double dt = atof(argv[3]);
    run(steps, dt);

	return 0;
}

