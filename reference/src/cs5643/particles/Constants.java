package cs5643.particles;

/**
 * Default constants. Add your own as necessary.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public class Constants {
	/** Mass of a particle. */
	public static double PARTICLE_MASS = 1.0;
	
	/** Radius of particle's sphere graphic. */
	public static double PARTICLE_RADIUS = 0.015;
	
	/** Camera rotation speed constants. */
	public static double CAM_SIN_THETA = Math.sin(0.2);
	public static double CAM_COS_THETA = Math.cos(0.2);
	
	/** Kernel Size */
	public static double KERNEL_SIZE = 0.1;
	
	/** Target Particle Density */
	public static double PARTICLE_DENSITY = 6378.0 * 100.0;
	
	/** Constraint Force Mixing Epsilon (to prevent divide by 0) */
	public static double CFM_EPSILON = 600;
	
	/** Artificial Pressure Strength */
	public static double ARTIFICIAL_PRESSURE_STRENGTH = 1e-5;
	
	/** Artificial Pressure Radius */
	public static double ARTIFICIAL_PRESSURE_RADIUS = 0.3 * KERNEL_SIZE;
	
	/** Artificial Pressure Power */
	public static double ARTIFICIAL_PRESSURE_POWER = 4;
	
	/** Artificial Viscosity */
	public static double ARTIFICIAL_VISCOSITY = 1e-2;
	
	/** Vorticity Constant */
	public static double VORTICITY_COEFFICIENT = 1e-4;
	
	/** Boundry Size */
	public static double BOUNDRY_SIZE = 1.0;
	
	static {
		//setHigh();
		setVeryHigh();
	}
	
	public static void setDefault() {
		PARTICLE_MASS = 1.0;
		PARTICLE_RADIUS = 0.015;
		KERNEL_SIZE = 0.1;
		PARTICLE_DENSITY = 6378.0 * 100.0;
		CFM_EPSILON = 600;
		ARTIFICIAL_PRESSURE_STRENGTH = 1e-5;
		ARTIFICIAL_PRESSURE_RADIUS = 0.3 * KERNEL_SIZE;
		ARTIFICIAL_PRESSURE_POWER = 4;
		ARTIFICIAL_VISCOSITY = 1e-2;
		VORTICITY_COEFFICIENT = 1e-4;
		BOUNDRY_SIZE = 1.0;
	}
	
	public static void setHigh() {
		PARTICLE_MASS = 1.0;
		PARTICLE_RADIUS = 0.005;
		KERNEL_SIZE = 0.1;
		PARTICLE_DENSITY = 6378.0 * 10000.0;
		CFM_EPSILON = 600;
		ARTIFICIAL_PRESSURE_STRENGTH = 1e-5;
		ARTIFICIAL_PRESSURE_RADIUS = 0.3 * KERNEL_SIZE;
		ARTIFICIAL_PRESSURE_POWER = 4;
		ARTIFICIAL_VISCOSITY = 1e-2;
		VORTICITY_COEFFICIENT = 1e-3;
		BOUNDRY_SIZE = 1.0;
	}
	
	public static void setVeryHigh() {
		PARTICLE_MASS = 1.0;
		PARTICLE_RADIUS = 0.005;
		KERNEL_SIZE = 0.1;
		PARTICLE_DENSITY = 6378.0 * 1000.0;
		CFM_EPSILON = 600;
		ARTIFICIAL_PRESSURE_STRENGTH = 1e-5;
		ARTIFICIAL_PRESSURE_RADIUS = 0.3 * KERNEL_SIZE;
		ARTIFICIAL_PRESSURE_POWER = 4;
		ARTIFICIAL_VISCOSITY = 1e-2;
		VORTICITY_COEFFICIENT = 1e-3;
		BOUNDRY_SIZE = 10.0;
	}
	
	/** Some flags */
	
	/** Prevent things from exploding due to compressed initial states */
	public static final boolean USE_INITIAL_VELOCITY_DAMPING = false;
	
	/** Use of scorr to simulate surface tension */
	public static final boolean USE_SURFACE_TENSION = true;
	
	/** Use vorticity confinement to maintain energy */
	public static final boolean USE_VORTICITY_CONFINEMENT = true;
	
	/** Put a roof on the box */
	public static final boolean USE_BOX_CEILING = false;
	
	/** Parallel computations. Set to 0 for serial. */
	public static int USE_CORES = 6;
}
