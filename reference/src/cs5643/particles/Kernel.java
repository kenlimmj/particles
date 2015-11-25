package cs5643.particles;

import egl.math.Vector3d;

public abstract class Kernel {
	//> Evaluate the kernel
	public abstract double evaluate(Vector3d r, double h);
	
	public double evaluate(Vector3d r) {
		return evaluate(r, Constants.KERNEL_SIZE);
	}
	
	public double evaluate(Vector3d p, Vector3d q) {
		return evaluate(p, q, Constants.KERNEL_SIZE);
	}
	
	public double evaluate(Vector3d p, Vector3d q, double h) {
		return evaluate(p.clone().sub(q), h);
	}
}
