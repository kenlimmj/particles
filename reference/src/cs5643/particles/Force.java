package cs5643.particles;

import javax.media.opengl.*;

/**
 * Particle system force.
 *
 * @author Doug James, January 2007
 */
public interface Force {
	/**
	 * Causes force to be applied to affected particles.
	 */
	public void applyForce(Particle p);
	
	/** Display any instructive force information, e.g., direction. */
	public void display(GL2 gl);
}
