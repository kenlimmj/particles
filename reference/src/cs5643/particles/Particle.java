package cs5643.particles;

import java.util.ArrayList;
import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import egl.math.Vector3d;

/**
 * Simple particle implementation, with miscellaneous adornments.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public class Particle {
	/** Display list index. */
	private static int PARTICLE_DISPLAY_LIST = -1;
	
	/** Highlighted appearance if true, otherwise white. */
	private boolean highlight = false;
	
	/** Default mass. */
	double m = Constants.PARTICLE_MASS;
	
	/** Deformed Position. */
	Vector3d x = new Vector3d();
	
	/** Undeformed/material Position. */
	Vector3d x0 = new Vector3d();
	
	/** Velocity. */
	Vector3d v = new Vector3d();
	
	/** Force accumulator. */
	Vector3d f = new Vector3d();
	
	/** Neighbors */
	ArrayList<Particle> neighbors = new ArrayList<Particle>();

	/** Candidate Position */
	Vector3d xc = new Vector3d();
	
	/** Candidate velocity */
	Vector3d vc = new Vector3d();
	
	/** Lambda */
	double lm = 0.0;
	
	/** Pressure */
	Vector3d dp = new Vector3d();
	
	/** Vorticity */
	Vector3d vo = new Vector3d();
	
	/** Flag for rigid particles */
	boolean rigid = false;
	
	/**
	 * Constructs particle with the specified material/undeformed
	 * coordinate, x0.
	 */
	Particle(Vector3d x0) {
		this.x0.set(x0);
		x.set(x0);
	}
	
	/** Draws spherical particle using a display list. */
	public void display(GL2 gl) {
		if (PARTICLE_DISPLAY_LIST < 0) {// MAKE DISPLAY LIST:
			int displayListIndex = gl.glGenLists(1);
			GLU glu = GLU.createGLU();
			GLUquadric quadric = glu.gluNewQuadric();
			gl.glNewList(displayListIndex, GL2.GL_COMPILE);
			glu.gluSphere(quadric, Constants.PARTICLE_RADIUS, 16, 8);
			gl.glEndList();
			glu.gluDeleteQuadric(quadric);
			
			// For older versions of JOGL
			//glu.destroy();
			
			System.out.println("MADE DISPLAY LIST " + displayListIndex + " : " + gl.glIsList(displayListIndex));
			PARTICLE_DISPLAY_LIST = displayListIndex;
		}
		
		/// COLOR: DEFAULT CYAN; GREEN IF HIGHLIGHTED
		float[] c = new float[] { 0f, 1f, 1f, 1f };//default: cyan
		if (highlight) {
			c[2] = 0;
		}
		
		// Hack to make things more colorful/interesting
		c[1] = (float) (x.y/Constants.BOUNDRY_SIZE);
		
		gl.glMaterialfv(GL2.GL_FRONT, GL2.GL_DIFFUSE, c, 0); // Color used by shader
		
		/// DRAW ORIGIN-CIRCLE TRANSLATED TO "p":
		gl.glPushMatrix();
		gl.glTranslated(x.x/Constants.BOUNDRY_SIZE, x.y/Constants.BOUNDRY_SIZE, x.z/Constants.BOUNDRY_SIZE);
		gl.glCallList(PARTICLE_DISPLAY_LIST); // Draw the particle
		gl.glPopMatrix();
	}
	
	/** Specifies whether particle should be drawn highlighted. */
	public void setHighlight(boolean highlight) {
		this.highlight = highlight;
	}
	
	/** True if particle should be drawn highlighted. */
	public boolean getHighlight() {
		return highlight;
	}
}
