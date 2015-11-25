package cs5643.particles;

import egl.math.Vector3d;

public abstract class KernelGradient {
	//> Evaluate the kernel gradient
	public abstract Vector3d evaluate(Vector3d r, double h);
	
	public Vector3d evaluate(Vector3d r) {
		return evaluate(r, Constants.KERNEL_SIZE);
	}
	
	public Vector3d evaluate(Vector3d p, Vector3d q) {
		return evaluate(p, q, Constants.KERNEL_SIZE);
	}
	
	public Vector3d evaluate(Vector3d p, Vector3d q, double h) {
		return evaluate(p.clone().sub(q), h);
	}
}
