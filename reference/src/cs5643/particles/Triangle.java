package cs5643.particles;

import javax.media.opengl.GL2;
import egl.math.Vector3d;

public class Triangle {
	
	Vector3d v0, u, v, n;
	
	public Triangle(Vector3d a, Vector3d b, Vector3d c) {
		v0 = a;
		u = new Vector3d(b).sub(a);
		v = new Vector3d(c).sub(a);
		n = new Vector3d(u).cross(v).normalize();
	}
	
	/** 
	 * Check for an intersection between this triangle and a line segment.
	 * Returns null if no intersection, or the normal if there is.
	 * 
	 * Adapted from http://geomalgorithms.com/a06-_intersect-2.html
	 */
	public Vector3d[] intersect(Vector3d lineStart, Vector3d lineEnd) {
		// Find segment - plane intersection
		Vector3d dir = lineEnd.clone().sub(lineStart);
		double r = n.dot(v0.clone().sub(lineStart)) / n.dot(dir);
		
		if (Double.isNaN(r) || r < 0 || 1 < r) {
			return null;
		}
		
		Vector3d p = dir.mul(r).add(lineStart);
		
		// Check if intersection is within triangle
		Vector3d w = new Vector3d(p).sub(v0);
		
		double uu = u.dot(u);
		double vv = v.dot(v);
		double uv = u.dot(v);
		double wu = w.dot(u);
		double wv = w.dot(v);
		
		double s = (uv * wv - vv * wu) / (uv * uv - uu * vv);
		double t = (uv * wu - uu * wv) / (uv * uv - uu * vv);
		
		if (s < 0 || t < 0 || 1 - s - t < 0) {
			return null;
		}
		
		return new Vector3d[] { p, n.clone() };
	}
	
	public void display(GL2 gl) {
		
	}
}
