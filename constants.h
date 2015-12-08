#define ITERATIONS_PER_STEP 10

#define GRAVITY -9.807
#define PI 3.14159265

#define BOX_SIZE 1.0
#define PARTICLE_DENSITY (6378.0 * 10000.0)
#define PARTICLE_MASS 1.0
#define FORCE_EPSILON 600.0
#define KERNEL_SIZE 0.1

#define ARTIFICIAL_PRESSURE_STRENGTH 1e-5
#define ARTIFICIAL_PRESSURE_RADIUS (0.3 * KERNEL_SIZE)
#define ARTIFICIAL_PRESSURE_POWER 4

#define ARTIFICIAL_VISCOSITY 1e-2
#define VORTICITY_COEFFICIENT 1e-3
#define USE_VORTICITY_CONFINEMENT 1
