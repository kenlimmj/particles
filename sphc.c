#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"

// Current time
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
  
    // IS THIS THE WRONG C Constant? Should be 7.5 / PI?
    const double c = 15.0 / PI;
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
    const double c = 315.0 / (64.0 * PI);
    double dx = cvx[i] - cvx[j];
    double dy = cvy[i] - cvy[j];
    double dz = cvz[i] - cvz[j];
    double r = sqrt(dx * dx + dy * dy + dz * dz);

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

    // TODO Find neighbors

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

            if (cpx[i] < -BOX_SIZE) { // Wall
                cpx[i] = -BOX_SIZE;
                cvx[i] = -cvx[i];
            }
            if (cpx[i] > BOX_SIZE) { // Wall
                cpx[i] = BOX_SIZE;
                cvx[i] = -cvx[i];
            }
            if (cpz[i] < -BOX_SIZE) { // Wall
                cpz[i] = -BOX_SIZE;
                cvz[i] = -cvz[i];
            }
            if (cpz[i] > BOX_SIZE) { // Wall
                cpz[i] = BOX_SIZE;
                cvz[i] = -cvz[i];
            }
        }
    }

    // TODO Vorticity confinement
    // Is the vorticity confinement the next two i to n loops?
    // Compute angular velocity
    for (int i = 0; i < n; ++i) {
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        double gv[3];
        double vx, vy, vz;
        for (int j = 0, l = nc[i]; j < l; ++j) {
            int nb = nbs[i][j];
          // I THINK VX, VY, VZ SHOULDNT BE VELOCITY BUT A DIFFERENT VARIABLE FOR VORTICITY
            vx = cvx[nb] - cvx[i];
            vy = cvy[nb] - cvy[i];
            vz = cvz[nb] - cvz[i];

            spikykernel(i, nb, gv);
          // SAME QUESTION HERE FOR VX VY VZ
            sx += vy * gv[2] - vz * gv[1];
            sy += vz * gv[0] - vx * gv[2];
            sz += vx * gv[1] - vy * gv[0];
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
        // DONT NEED VX, VY, VZ??
        double vx, vy, vz;
        double l;
        for (int j = 0, l = nc[i]; j < l; ++j) {
            int nb = nbs[i][j];
            spikykernel(i, nb, gv);
            // USED L TWICE IN FOR LOOP AND IN BELOW
            l = sqrt(dpx[nb] * dpx[nb] + dpy[nb] * dpy[nb] + dpz[nb] * dpz[nb]);
            gv[0] *= l / PARTICLE_DENSITY;
            gv[1] *= l / PARTICLE_DENSITY;
            gv[2] *= l / PARTICLE_DENSITY;

            sx += gv[0];
            sy += gv[1];
            sz += gv[2];
        }
        l = sqrt(sx * sx + sy * sy + sz * sz);
        if (l < 1e-5) {
            l = 1e-5;
        }
        sx /= l;
        sy /= l;
        sz /= l;

        vox[i] = (sy * dpz[i] - sz * dpy[i]) * VORTICITY_COEFFICIENT;
        voy[i] = (sz * dpx[i] - sx * dpz[i]) * VORTICITY_COEFFICIENT;
        // I THINK BELOW SHOULD BE SX*DPY - SY*DPX
        voz[i] = (sx * dpy[i] - sy * dpz[i]) * VORTICITY_COEFFICIENT;
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
    dt = 0.0;
    for (int i = 0; i < 10; ++i) {
        step();
        time += dt;
    }
}

void init() {

}

int main(int argc, char ** argv) {
	return 0;
}
