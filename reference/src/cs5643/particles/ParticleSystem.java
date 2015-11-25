package cs5643.particles;

import java.awt.Font;
import java.util.*;
import javax.media.opengl.*;
import com.jogamp.opengl.util.awt.TextRenderer;
import com.jogamp.opengl.util.glsl.*;
import egl.math.Vector3d;

/**
 * Maintains dynamic lists of Particle and Force objects, and provides
 * access to their state for numerical integration of dynamics.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public class ParticleSystem {
	/** Current simulation time. */
	public double time = 0;
	
	/** List of Particle objects. */
	public ArrayList<Particle> P = new ArrayList<Particle>();
	
	/** List of rigid triangles */
	public ArrayList<Triangle> T = new ArrayList<Triangle>();
	
	/** List of Force objects. */
	public ArrayList<Force> F = new ArrayList<Force>();
	
	/**
	 * true iff prog has been initialized. This cannot be done in the
	 * constructor because it requires a GL2 reference.
	 */
	private boolean init = false;
	
	/** Filename of vertex shader source. */
	public static final String[] VERT_SOURCE = { "vert.glsl" };
	
	/** Filename of fragment shader source. */
	public static final String[] FRAG_SOURCE = { "frag.glsl" };
	
	/** The shader program used by the particles. */
	ShaderProgram prog;
	
	/** Number of iterations per step */
	public static int ITERATIONS_PER_STEP = 10;
	
	/** Kernel */
	Kernel poly6, spiky, viscosity;
	KernelGradient spikyGrad;
	
	private int DEBUG = 0;
	
	/** Cores */
	public static int CORES = 4;
	
	/** Basic constructor. */
	public ParticleSystem() {
		// Initialize kernels
		poly6 = new Kernel() {
			double c = 315.0 / (64.0 * Math.PI);
			
			@Override
			public double evaluate(Vector3d v, double h) {
				double r = v.len();
				
				return (r > h) ? 0 : c * Math.pow(h * h - r * r, 3) / Math.pow(h, 9);
			}
		};
		
		spiky = new Kernel() {
			double c = 15.0 / Math.PI;
			
			@Override
			public double evaluate(Vector3d v, double h) {
				double r = v.len();
				
				return (r > h) ? 0 : c * Math.pow(h - r, 3) / Math.pow(h, 6);
			}
		};
		
		spikyGrad = new KernelGradient() {
			double c = 3.0 * 15.0 / Math.PI;
			
			@Override
			public Vector3d evaluate(Vector3d v, double h) {
				double r = v.len();
				
				if (r < 1e-4) r = 1e-4;
				
				double f = c / Math.pow(h, 6) * Math.pow(h - r, 2) / r;
				
				return (r > h) ? new Vector3d() : new Vector3d(f).mul(v);
			}
		};
		
		viscosity = new Kernel() {
			double c = 15.0 / (2 * Math.PI);
			
			@Override
			public double evaluate(Vector3d v, double h) {
				double r = v.len();
				
				if (r < 1e-4) r = 1e-4;
				
				return (r > h) ? 0 : c * Math.pow(h, 3)
					* (-(r * r * r) / (2 * h * h * h) + (r * r) / (h * h) + (h) / (2 * r) - 1);
			}
		};
		
		PRESSURE_RADIUS_FACTOR = 1.0 / poly6.evaluate(new Vector3d(0, 0, Constants.ARTIFICIAL_PRESSURE_RADIUS));
	}
	
	/**
	 * Set up the GLSL program. This requires that the current directory (i.e. the package in which
	 * this class resides) has a vertex and fragment shader.
	 */
	public synchronized void init(GL2 gl) {
		if (init) return;
		
		prog = new ShaderProgram();
		ShaderCode vert_code = ShaderCode.create(gl, GL2ES2.GL_VERTEX_SHADER, 1, this.getClass(), VERT_SOURCE, false);
		ShaderCode frag_code = ShaderCode.create(gl, GL2ES2.GL_FRAGMENT_SHADER, 1, this.getClass(), FRAG_SOURCE, false);
		
		System.out.println(vert_code);
		System.out.println(frag_code);
		
		if (!prog.add(gl, vert_code, System.err) || !prog.add(gl, frag_code, System.err)) {
			System.err.println("WARNING: shader did not compile");
			prog.init(gl); // Initialize empty program
		}
		else {
			prog.link(gl, System.err);
		}
		
		init = true;
	}
	
	/** Adds a force object (until removed) */
	public synchronized void addForce(Force f) {
		F.add(f);
	}
	
	/** Useful for removing temporary forces, such as user-interaction
	 * spring forces. */
	public synchronized void removeForce(Force f) {
		F.remove(f);
	}
	
	/** Creates particle and adds it to the particle system.
	 * @param p0 Undeformed/material position.
	 * @return Reference to new Particle.
	 */
	public synchronized Particle createParticle(Vector3d p0) {
		Particle newP = new Particle(p0);
		P.add(newP);
		return newP;
	}
	
	/**
	 * Helper-function that computes nearest particle to the specified
	 * (deformed) position.
	 * @return Nearest particle, or null if no particles.
	 */
	public synchronized Particle getNearestParticle(Vector3d x) {
		Particle minP = null;
		double minDistSq = Double.MAX_VALUE;
		for (Particle particle : P) {
			double distSq = x.distSq(particle.x);
			if (distSq < minDistSq) {
				minDistSq = distSq;
				minP = particle;
			}
		}
		return minP;
	}
	
	/** Create rigid particle */
	public synchronized void addTriangle(Triangle t) {
		T.add(t);
	}
	
	/** Moves all particles to undeformed/materials positions, and
	 * sets all velocities to zero. Synchronized to avoid problems
	 * with simultaneous calls to advanceTime(). */
	public synchronized void reset() {
		for (Particle p : P) {
			p.x.set(p.x0);
			p.v.set(0, 0, 0);
			p.f.set(0, 0, 0);
			p.setHighlight(false);
		}
		time = 0;
	}
	
	/** Hash three integers for the grid */
	public Integer hash(int x, int y, int z) {
		//if (x < 0 || x > 10 || y < 0 || y > 10 || z < 0 || z > 10) return -1;
		x = x + 100;
		y = y + 100;
		z = z + 100;
		return x + 400 * y + 40000 * z + 1;
	}
	
	/**
	 * Incomplete/Debugging implementation of Forward-Euler step.
	 */
	public synchronized void advanceTime(double dt) {
		long startTime = System.currentTimeMillis();
		
		System.out.print("Start");
		
		// Clear force accumulators
		clearForces();
		
		// Apply forces
		applyForces(F);
		
		// Compute candidate velocity and positions
		computeCandidates(dt);
		
		// Grid hashing
		Map<Integer,ArrayList<Particle>> grid = hashParticles();
		
		// Find neighbors using the kernel
		//*
		findNeighborsGrid(grid);
		/*/
		findNeighborsNaive();
		//*/
		
		System.out.print(" - Iterations");
		
		// Perform adjustment iterations
		int iteration = 0;
		while (iteration++ < ITERATIONS_PER_STEP) {
			// Compute lambda
			computeLambda();
			
			// Jiggle particles
			solveParticles(dt);
		}
		
		System.out.print(" - Vorticity");
		
		// Vorticity Confinement
		if (Constants.USE_VORTICITY_CONFINEMENT) {
			computeVorticityOmega();
			
			computeVorticityForce();
		}
		
		System.out.print(" - XSPH");
		
		// XSPH Viscosity
		computeXSPH();
		
		// Update position and velocity
		updateParticle();
		
		System.out.print(" - Done");
		
		long endTime = System.currentTimeMillis();
		System.out.println(" - " + (endTime - startTime) + " ms");
		
		time += dt;
	}
	
	/** Parallel lambda computation */
	private class Worker_Lambda implements Runnable {
		private int start, end;
		public Worker_Lambda(int start, int end) {
			this.start = start;
			if (end>P.size()) end = P.size();
			this.end = end;
		}
		
		public void run() {
			for(int i=start; i<end; i++) {
				Particle p = P.get(i);
				
				// Compute density
				double rho = 0;
				double numerator = 0;
				double denominator = 0;
				
				// Compute density
				for (Particle q : p.neighbors) {
					rho += q.m * poly6.evaluate(p.xc, q.xc);
				}
				numerator = rho / Constants.PARTICLE_DENSITY - 1;
				
				// Compute gradient
				Vector3d sum = new Vector3d(), grad;
				for (Particle q : p.neighbors) {
					// k = i
					grad = spikyGrad.evaluate(p.xc, q.xc);
					sum.add(grad);
					
					// k = j
					grad.div(Constants.PARTICLE_DENSITY).negate();
					denominator += grad.lenSq();
				}
				sum.div(Constants.PARTICLE_DENSITY);
				denominator += sum.lenSq();
				
				p.lm = -numerator / (denominator + Constants.CFM_EPSILON);
			}
		}
	}
	
	/** Parallel particle solver */
	private class Worker_Solve implements Runnable {
		private int start, end;
		private double dt;
		public Worker_Solve(int start, int end, double dt) {
			this.start = start;
			if (end>P.size()) end = P.size();
			this.end = end;
			this.dt = dt;
		}
		
		public void run() {
			for(int i=start; i<end; i++) {
				Particle p = P.get(i);
				
				computePressure(p);
				
				// Update candidate position and velocity
				updateCandidate(p, dt);
				
				// Collisions
				computeCollision(p);
			}
		}
	}
	
	/** ===== PBF Parts ===== */
	
	/** Clear accumulated forces on the particle */
	public void clearForces() {
		for (Particle p : P) {
			p.f.set(0, 0, 0);
		}
	}
	
	/** Apply force to this particle system */
	public void applyForces(ArrayList<Force> F) {
		for (Particle p : P) {
			for (Force f : F) {
				f.applyForce(p);
			}
			
			// Vorticity force
			p.f.add(p.vo);
		}
	}
	
	/** Compute candidate velocity and position */
	public void computeCandidates(double dt) {
		for (Particle p : P) {
			p.vc.set(p.v).addMultiple(dt, p.f);
			if (Constants.USE_INITIAL_VELOCITY_DAMPING && time < 1.0 && p.vc.len() > 2) {
				p.vc.normalize();
			}
			p.xc.set(p.x).addMultiple(dt, p.vc);
		}
	}
	
	/** Hash the particles into a grid */
	public Map<Integer,ArrayList<Particle>> hashParticles() {
		Map<Integer,ArrayList<Particle>> grid = new HashMap<Integer,ArrayList<Particle>>(150);
		for (Particle p : P) {
			int x = (int) (p.xc.x / Constants.KERNEL_SIZE);
			int y = (int) (p.xc.y / Constants.KERNEL_SIZE);
			int z = (int) (p.xc.z / Constants.KERNEL_SIZE);
			int hash = hash(x, y, z);
			if (grid.containsKey(hash)) {
				grid.get(hash).add(p);
			}
			else {
				ArrayList<Particle> list = new ArrayList<Particle>();
				list.add(p);
				grid.put(hash, list);
			}
		}
		
		return grid;
	}
	
	/** Possible grid directions for neighbors */
	final int[][] directions = {
		{ -1, -1, -1 }, { -1, -1, 0 }, { -1, -1, 1 }, { -1, 0, -1 }, { -1, 0, 0 }, { -1, 0, 1 }, { -1, 1, -1 },
		{ -1, 1, 0 }, { -1, 1, 1 },
		{ 0, -1, -1 }, { 0, -1, 0 }, { 0, -1, 1 }, { 0, 0, -1 }, { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, -1 }, { 0, 1, 0 },
		{ 0, 1, 1 },
		{ 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 }, { 1, 0, -1 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 1, -1 }, { 1, 1, 0 },
		{ 1, 1, 1 }
	};
	
	/** Find neighbors using the hashed grid */
	public void findNeighborsGrid(Map<Integer,ArrayList<Particle>> grid) {
		for (Particle p : P) {
			p.neighbors.clear();
			int x = (int) (p.xc.x / Constants.KERNEL_SIZE);
			int y = (int) (p.xc.y / Constants.KERNEL_SIZE);
			int z = (int) (p.xc.z / Constants.KERNEL_SIZE);
			
			for (int i = 0; i < directions.length; i++) { // Nearby grid spaces
				int hash = hash(x + directions[i][0], y + directions[i][1], z + directions[i][2]);
				if (grid.containsKey(hash)) {
					for (Particle q : grid.get(hash)) {
						if (p.xc.dist(q.xc) <= Constants.KERNEL_SIZE) {
							p.neighbors.add(q);
						}
					}
				}
			}
		}
	}
	
	/** Find neighbors using the naive N^2 method */
	public void findNeighborsNaive() {
		for (Particle p : P) {
			p.neighbors.clear();
			for (Particle q : P) {
				if (p.xc.dist(q.xc) <= Constants.KERNEL_SIZE) {
					p.neighbors.add(q);
				}
			}
		}
	}
	
	/** Compute the density of each particle and the multiplier lambda */
	public void computeLambda() {
		if (Constants.USE_CORES>0) {
			int workers = Constants.USE_CORES;
			int workSize = P.size()/workers + 1;
			int offset = 0;
			
			Thread[] threads = new Thread[workers];
			for (int i=0; i<workers; i++) {
				Worker_Lambda w = new Worker_Lambda(offset, offset + workSize);
				threads[i] = new Thread(w);
				threads[i].start();
				offset += workSize;
			}
			
			for (int i=0; i<workers; i++) {
				try {
					threads[i].join();
				}
				catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		else {
			for (Particle p : P) {
				// Compute density
				double rho = 0;
				double numerator = 0;
				double denominator = 0;
				
				// Compute density
				for (Particle q : p.neighbors) {
					rho += q.m * poly6.evaluate(p.xc, q.xc);
				}
				numerator = rho / Constants.PARTICLE_DENSITY - 1;
				
				// Compute gradient
				Vector3d sum = new Vector3d(), grad;
				for (Particle q : p.neighbors) {
					// k = i
					grad = spikyGrad.evaluate(p.xc, q.xc);
					sum.add(grad);
					
					// k = j
					grad.div(Constants.PARTICLE_DENSITY).negate();
					denominator += grad.lenSq();
				}
				sum.div(Constants.PARTICLE_DENSITY);
				denominator += sum.lenSq();
				
				p.lm = -numerator / (denominator + Constants.CFM_EPSILON);
			}
		}
	}
	
	/** Solve particle constraints */
	public void solveParticles(double dt) {
		if (Constants.USE_CORES>0) {
			int workers = Constants.USE_CORES;
			int workSize = P.size()/workers + 1;
			int offset = 0;
			
			Thread[] threads = new Thread[workers];
			for (int i=0; i<workers; i++) {
				Worker_Solve w = new Worker_Solve(offset, offset + workSize, dt);
				threads[i] = new Thread(w);
				threads[i].start();
				offset += workSize;
			}
			
			for (int i=0; i<workers; i++) {
				try {
					threads[i].join();
				}
				catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		else {
			for (Particle p : P) {
				computePressure(p);
				
				// Update candidate position and velocity
				updateCandidate(p, dt);
				
				// Collisions
				computeCollision(p);
			}
		}
	}
	
	public static double PRESSURE_RADIUS_FACTOR;
	
	/** Compute the change in position due to pressure */
	public void computePressure(Particle p) {
		// Particle pressures
		Vector3d sum = new Vector3d();
		for (Particle q : p.neighbors) {
			double scorr = 0;
			if (Constants.USE_SURFACE_TENSION) {
				scorr = -Constants.ARTIFICIAL_PRESSURE_STRENGTH
					* Math.pow(
						poly6.evaluate(p.xc, q.xc) * PRESSURE_RADIUS_FACTOR,
						Constants.ARTIFICIAL_PRESSURE_POWER
						);
			}
			
			sum.add(
				spikyGrad.evaluate(p.xc, q.xc)
					.mul(p.lm + q.lm + scorr)
				);
			if (DEBUG > 0) System.out.println(spikyGrad.evaluate(p.xc, q.xc));
		}
		p.dp.set(sum).div(Constants.PARTICLE_DENSITY);
	}
	
	/** Update candidate positions due to pressure */
	public void updateCandidate(Particle p, double dt) {
		p.xc.add(p.dp);
		p.vc.set(p.xc).sub(p.x).div(dt);
	}
	
	/** Handle collisions */
	public void computeCollision(Particle p) {
		final double size = Constants.BOUNDRY_SIZE;
		// Box
		if (p.xc.y < 0) { // Floor
			p.xc.y = 0;
			p.vc.y = -p.vc.y;
		}
		if (Constants.USE_BOX_CEILING) { // Ceiling
			if (p.xc.y > size) {
				p.xc.y = size;
				p.vc.y = 0;
			}
		}
		if (p.xc.x < 0) { // Walls
			p.xc.x = 0;
			p.vc.x = -p.vc.x;
		}
		if (p.xc.x > size) { // Walls
			p.xc.x = size;
			p.vc.x = -p.vc.x;
		}
		if (p.xc.z < 0) { // Walls
			p.xc.z = 0;
			p.vc.z = -p.vc.z;
		}
		if (p.xc.z > size) { // Walls
			p.xc.z = size;
			p.vc.z = -p.vc.z;
		}
		
		// Triangles
		Vector3d[] intersection;
		for (Triangle t : T) {
			intersection = t.intersect(p.x, p.xc);
			if (intersection != null) {
				p.xc.set(intersection[0]);
				p.vc.sub(intersection[1].mul(p.vc.dot(intersection[1])));
			}
		}
	}
	
	/** Compute vorticity omega */
	public void computeVorticityOmega() {
		for (Particle p : P) {
			// Compute vorticity
			Vector3d sum = new Vector3d();
			for (Particle q : p.neighbors) {
				Vector3d vel = new Vector3d(q.vc).sub(p.vc);
				Vector3d grad = spikyGrad.evaluate(p.xc, q.xc);
				sum.add(vel.cross(grad));
			}
			// Reusing the variable dp just because
			p.dp.set(sum);
		}
	}
	
	/** Compute vorticity force */
	public void computeVorticityForce() {
		for (Particle p : P) {
			// Compute vorticity forces
			Vector3d sum = new Vector3d();
			for (Particle q : p.neighbors) {
				sum.add(
					spikyGrad.evaluate(p.xc, q.xc)
						.mul(q.dp.len()).div(Constants.PARTICLE_DENSITY));
			}
			
			Vector3d N = normalize(sum);
			Vector3d vf = N.cross(p.dp).mul(Constants.VORTICITY_COEFFICIENT);
			p.vo.set(vf);
		}
	}
	
	/** Normalize a vector and handle small case */
	private Vector3d normalize(Vector3d v) {
		double l = v.len();
		if (l < 1e-5) {
			l = 1e-5;
		}
		return new Vector3d(v).div(l);
	}
	
	/** Compute XSPH velocities */
	public void computeXSPH() {
		for (Particle p : P) {
			Vector3d weightedSum = new Vector3d();
			for (Particle q : p.neighbors) {
				Vector3d diff = new Vector3d(q.vc).sub(p.vc);
				weightedSum.addMultiple(poly6.evaluate(p.vc, q.vc), diff);
			}
			p.vc.addMultiple(Constants.ARTIFICIAL_VISCOSITY, weightedSum);
		}
	}
	
	/** Finalize all updates to particles */
	public void updateParticle() {
		for (Particle p : P) {
			p.v.set(p.vc);
			p.x.set(p.xc);
		}
	}
	
	/** Text renderer object */
	TextRenderer textrenderer;
	
	/**
	 * Displays Particle and Force objects.
	 */
	public synchronized void display(GL2 gl) {
		for (Force force : F) {
			force.display(gl);
		}
		
		if (!init) {
			textrenderer = new TextRenderer(new Font("SansSerif", Font.PLAIN, 16));
			init(gl);
		}
		
		prog.useProgram(gl, true);
		
		for (Particle p : P) {
			p.display(gl);
		}
		
		for (Triangle t : T) {
			t.display(gl);
		}
		
		prog.useProgram(gl, false);
		
		textrenderer.beginRendering(800, 600, true);
		textrenderer.draw(String.valueOf(time), 10, 10);
		textrenderer.endRendering();
	}
}
