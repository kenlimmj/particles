#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "constants.h"

#define ITERATIONS_PER_STEP 10
#define USE_VORTICITY_CONFINEMENT 1

void clearForces() {

}

void applyForces(double[] forces) {

}

void computeCandidates(double dt) {

}

grid hashParticles() {

}

void findNeighbors(grid grid) {

}

void computeLambda() {

}

void solveParticles(double dt) {

}

void computeVorticityOmega() {

}

void computeVorticityForce() {

}

void computeXSPH() {

}

void updateParticle() {

}

void advanceTime(double dt) {
  double startTime = omp_get_time();

  clearForces();
  applyForces(F);

  computeCandidates(dt);

  grid = hashParticles();

  findNeighbors(grid);

  int i = 0;
  while (i < ITERATIONS_PER_STEP) {
    computeLambda();
    solveParticles(dt);
  }

  if (USE_VORTICITY_CONFINEMENT) {
    computeVorticityOmega();
    computeVorticityForce();
  }

  computeXSPH();
  updateParticle();

  double endTime = omp_get_time();

  time += dt;
}

int main() {
  
  return 0;
}
