package cs5643.particles;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.io.File;
import java.io.IOException;
import java.util.Random;
import javax.swing.*;
import javax.imageio.ImageIO;
import javax.media.opengl.*;
import javax.media.opengl.awt.GLCanvas;
import javax.media.opengl.glu.*;
import com.jogamp.opengl.util.*;
import egl.math.Vector2d;
import egl.math.Vector3d;

/**
 * CS5643: Assignment #1: Smoothed-Particle Hydrodynamics
 *
 * main() entry point class that initializes ParticleSystem, OpenGL
 * rendering, and GUI that manages GUI/mouse events.
 *
 * Spacebar toggles simulation advance.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public class ParticleSystemBuilder implements GLEventListener {
	private FrameExporter frameExporter;
	
	private static int N_STEPS_PER_FRAME = 5;
	
	private GLU glu;
	
	/** Default graphics time step size. */
	public static final double DT = 0.01;
	
	/** Main window frame. */
	JFrame frame = null;
	
	private int width, height;
	
	/** The single ParticleSystem reference. */
	ParticleSystem PS;
	
	/** Object that handles all GUI and user interactions of building
	 * Task objects, and simulation. */
	BuilderGUI gui;
	
	/** Position of the camera. */
	public Vector3d eyePos = new Vector3d(14, 10, 10);
	
	/** Position of the camera's focus. */
	public Vector3d targetPos = new Vector3d(0.5, 0.5, 0.5);
	
	/** Position of the light. Fixed at the location of the camera. */
	private float[] lightPos = { 0f, 0f, 0f, 1f };
	
	/** Main constructor. Call start() to begin simulation. */
	ParticleSystemBuilder() {
		PS = new ParticleSystem();
	}
	
	/**
	 * Builds and shows windows/GUI, and starts simulator.
	 */
	public void start() {
		if (frame != null) return;
		
		gui = new BuilderGUI();
		
		frame = new JFrame("CS567 Particle System Builder");
		GLProfile glp = GLProfile.getDefault();
		GLCapabilities glc = new GLCapabilities(glp);
		GLCanvas canvas = new GLCanvas(glc);
		canvas.addGLEventListener(this);
		frame.add(canvas);
		
		canvas.addMouseListener(gui);
		canvas.addMouseMotionListener(gui);
		canvas.addKeyListener(gui);
		
		final Animator animator = new Animator(canvas);
		frame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				// Run this on another thread than the AWT event queue to
				// make sure the call to Animator.stop() completes before
				// exiting
				new Thread(new Runnable() {
					public void run() {
						animator.stop();
						System.exit(0);
					}
				}).start();
			}
		});
		
		frame.pack();
		frame.setSize(818, 647);
		frame.setLocation(360, 0);
		frame.setVisible(true);
		animator.start();
	}
	
	/** GLEventListener implementation: Initializes JOGL renderer. */
	public void init(GLAutoDrawable drawable) {
		// DEBUG PIPELINE (can use to provide GL error feedback... disable for speed)
		//drawable.setGL(new DebugGL(drawable.getGL()));
		
		GL2 gl = drawable.getGL().getGL2();
		System.err.println("INIT GL IS: " + gl.getClass().getName());
		
		gl.setSwapInterval(1);
		gl.glEnable(GL2.GL_DEPTH_TEST);
		gl.glLineWidth(1);
		
		gl.glEnable(GL2.GL_NORMALIZE);
		
		// SETUP LIGHTING
		float[] lightAmbient = { 0f, 0f, 0f, 1f };
		float[] lightDiffuse = { 0.9f, 0.9f, 0.9f, 1f };
		float[] lightSpecular = { 1f, 1f, 1f, 1f };
		
		gl.glMatrixMode(GL2.GL_MODELVIEW);
		gl.glLoadIdentity();
		gl.glLightfv(GL2.GL_LIGHT0, GL2.GL_POSITION, lightPos, 0);
		gl.glLightfv(GL2.GL_LIGHT0, GL2.GL_AMBIENT, lightAmbient, 0);
		gl.glLightfv(GL2.GL_LIGHT0, GL2.GL_DIFFUSE, lightDiffuse, 0);
		gl.glLightfv(GL2.GL_LIGHT0, GL2.GL_SPECULAR, lightSpecular, 0);
		gl.glEnable(GL2.GL_LIGHT0);
		
		gui.initialize();
	}
	
	/** GLEventListener implementation */
	public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged) {
	}
	
	/** GLEventListener implementation */
	public void dispose(GLAutoDrawable drawable) {
	}
	
	/** GLEventListener implementation */
	public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {
		System.out.println("width=" + width + ", height=" + height);
		height = Math.max(height, 1); // avoid height=0;
		
		this.width = width;
		this.height = height;
		
		GL2 gl = drawable.getGL().getGL2();
		gl.glViewport(0, 0, width, height);
		
	}
	
	/**
	 * Main event loop: OpenGL display + simulation
	 * advance. GLEventListener implementation.
	 */
	public void display(GLAutoDrawable drawable) {
		GL2 gl = drawable.getGL().getGL2();
		gl.glClearColor(0.1f, 0.1f, 0.2f, 1f);
		gl.glClear(GL2.GL_COLOR_BUFFER_BIT | GL2.GL_DEPTH_BUFFER_BIT);
		
		/// GET READY TO DRAW:
		gl.glMatrixMode(GL2.GL_MODELVIEW);
		if (glu == null) glu = GLU.createGLU();
		gl.glLoadIdentity();
		glu.gluPerspective(5, (float) width / height, 1, 100);
		glu.gluLookAt(eyePos.x, eyePos.y, eyePos.z, targetPos.x, targetPos.y, targetPos.z, 0, 1, 0);
		
		/// DRAW COMPUTATIONAL CELL BOUNDARY:
		gl.glColor3f(1, 0, 0);
		gl.glBegin(GL2.GL_LINE_LOOP);
		gl.glVertex3d(0, 0, 0);
		gl.glVertex3d(1, 0, 0);
		gl.glVertex3d(1, 1, 0);
		gl.glVertex3d(0, 1, 0);
		gl.glEnd();
		gl.glBegin(GL2.GL_LINE_LOOP);
		gl.glVertex3d(0, 0, 1);
		gl.glVertex3d(1, 0, 1);
		gl.glVertex3d(1, 1, 1);
		gl.glVertex3d(0, 1, 1);
		gl.glEnd();
		gl.glBegin(GL2.GL_LINES);
		gl.glVertex3d(0, 0, 0);
		gl.glVertex3d(0, 0, 1);
		gl.glEnd();
		gl.glBegin(GL2.GL_LINES);
		gl.glVertex3d(1, 0, 0);
		gl.glVertex3d(1, 0, 1);
		gl.glEnd();
		gl.glBegin(GL2.GL_LINES);
		gl.glVertex3d(1, 1, 0);
		gl.glVertex3d(1, 1, 1);
		gl.glEnd();
		gl.glBegin(GL2.GL_LINES);
		gl.glVertex3d(0, 1, 0);
		gl.glVertex3d(0, 1, 1);
		gl.glEnd();
		
		/// SIMULATE/DISPLAY HERE (Handled by BuilderGUI):
		gui.simulateAndDisplayScene(gl);
	}
	
	/** Interaction central: Handles windowing/mouse events, and building state. */
	class BuilderGUI implements MouseListener, MouseMotionListener, KeyListener {
		boolean simulate = false;
		boolean single = false;
		
		/** Current build task (or null) */
		Task task;
		
		JFrame guiFrame;
		TaskSelector taskSelector = new TaskSelector();
		
		BuilderGUI() {
			guiFrame = new JFrame("Tasks");
			guiFrame.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
			guiFrame.setLayout(new SpringLayout());
			guiFrame.setLayout(new GridLayout(6, 1));
			
			/* Add new task buttons here, then add their functionality below. */
			ButtonGroup buttonGroup = new ButtonGroup();
			AbstractButton[] buttons = { new JButton("Reset"),
				new JButton("Load File"),
				new JToggleButton("Create Particle", false),
				new JButton("Toggle Gravity"),
				new JButton("Toggle Wind"),
			};
			
			for (int i = 0; i < buttons.length; i++) {
				buttonGroup.add(buttons[i]);
				guiFrame.add(buttons[i]);
				buttons[i].addActionListener(taskSelector);
			}
			
			guiFrame.setSize(320, 320);
			guiFrame.setVisible(true);
			
			task = null; // Set default task here
		}
		
		private Force gravityForce, windForce;
		
		/**
		 * Initialize some simulation variables.
		 */
		void initialize() {
			gravityForce = new Force() {
				@Override
				public void applyForce(Particle p) {
					//> Gravity is not affected by mass?
					p.f.sub(new Vector3d(0.0, 10.0, 0.0));
				}
				
				@Override
				public void display(GL2 gl) {
					// TODO Auto-generated method stub
					
				}
			};
			windForce = new Force() {
				Vector3d windVector = new Vector3d(4, 0, -4);
				
				@Override
				public void applyForce(Particle p) {
					Vector3d diff = new Vector3d(windVector).sub(p.v);
					p.f.add(diff);
				}
				
				@Override
				public void display(GL2 gl) {
					// TODO Auto-generated method stub
					
				}
			};
			
			PS.addForce(gravityForce);
			
			int c = 0;
			/*// Pineapple
			try {
				BufferedImage i = ImageIO.read(new File("Pineapple@64.bmp"));
				for (int x = 0; x < 64; x++) {
					for (int y = 0; y < 64; y++) {
						int col = i.getRGB(x, y);
						if (col < -1) {
							double dx = x * 0.06 + 0.6;
							double dy = -y * 0.06 + 5.0;
							double dz = 0.0;
							for (int z = 0; z < 20; z++) {
								dz = z * 0.06 + 1.9;
								PS.createParticle(new Vector3d(dx, dy, dz));
								c++;
							}
						}
					}
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			//*/
			
			//*// Pineapple x2
			try {
				BufferedImage i = ImageIO.read(new File("Pineapple@256.bmp"));
				for (int x = 0; x < 256; x++) {
					for (int y = 0; y < 256; y++) {
						int col = i.getRGB(x, y);
						if (col < -1) {
							double dx = x * 0.06 - 2.18 - 10.0;
							double dy = -y * 0.06 + 15.0;
							double dz = 0.0;
							for (int z = 0; z < 10; z++) {
								dz = z * 0.06 + 4.5;
								PS.createParticle(new Vector3d(dz, dy, -dx));
								c++;
							}
						}
					}
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			//*/
			
			/*// Cube
			double d = 0.06;
			for (double x = 0.1; x <= 2.1; x += d) {
				for (double y = 2.9; y <= 4.9; y += d) {
					for (double z = 0.1; z <= 2.1; z += d) {
						PS.createParticle(new Vector3d(x, y, z));
						c++;
					}
				}
			}
			//*/
			
			/*// Random
			Random r = new Random(1);
			for (int i=0; i<4000; i++) {
				PS.createParticle(new Vector3d(r.nextDouble(), r.nextDouble() * 0.6, r.nextDouble()));
				c++;
			}
			//*/
			
			System.out.println("Created " + c + " particles");
		}
		
		/** Simulate then display particle system and any builder
		 * adornments. */
		void simulateAndDisplayScene(GL2 gl) {
			if (simulate || single) {
				int nSteps = N_STEPS_PER_FRAME;
				double dt = DT / (double) nSteps;
				for (int k = 0; k < nSteps; k++) {
					PS.advanceTime(dt);
				}
				
				if (frameExporter != null) {
					frameExporter.writeFrame(gl);
				}
				
				single = false;
			}
			
			// Draw particles, forces, etc.
			PS.display(gl);
			
			// Display task if any
			if (task != null) task.display(gl);
		}
		
		/**
		 * ActionListener implementation to manage Task selection
		 * using (radio) buttons.
		 */
		class TaskSelector implements ActionListener {
			/**
			 * Resets ParticleSystem to undeformed/material state,
			 * disables the simulation, and removes the active Task.
			 */
			void resetToRest() {
				PS.reset(); //synchronized
				simulate = false;
				task = null;
			}
			
			/** Creates new Task objects to handle specified button action.
			 *  Switch to a new task, or perform custom button actions here.
			 */
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();
				System.out.println(cmd);
				
				if (cmd.equals("Reset")) {
					if (task != null) {
						task.reset();
					}
					else {
						resetToRest(); // set task=null
					}
				}
				else if (cmd.equals("Create Particle")) {
					task = new CreateParticleTask();
				}
				else if (cmd.equals("Load File")) {
					loadFrameFromFile();
				}
				else if (cmd.equals("Toggle Gravity")) {
					if (PS.F.contains(gravityForce)) {
						PS.removeForce(gravityForce);
					}
					else {
						PS.addForce(gravityForce);
					}
				}
				else if (cmd.equals("Toggle Wind")) {
					if (PS.F.contains(windForce)) {
						PS.removeForce(windForce);
					}
					else {
						PS.addForce(windForce);
					}
				}
				else {
					System.out.println("UNHANDLED ActionEvent: " + e);
				}
			}
		}
		
		// Methods required for the implementation of MouseListener
		public void mouseEntered(MouseEvent e) {
			if (task != null) task.mouseEntered(e);
		}
		
		public void mouseExited(MouseEvent e) {
			if (task != null) task.mouseExited(e);
		}
		
		public void mousePressed(MouseEvent e) {
			if (task != null) task.mousePressed(e);
		}
		
		public void mouseReleased(MouseEvent e) {
			if (task != null) task.mouseReleased(e);
		}
		
		public void mouseClicked(MouseEvent e) {
			if (task != null) task.mouseClicked(e);
		}
		
		// Methods required for the implementation of MouseMotionListener
		public void mouseDragged(MouseEvent e) {
			if (task != null) task.mouseDragged(e);
		}
		
		public void mouseMoved(MouseEvent e) {
			if (task != null) task.mouseMoved(e);
		}
		
		// Methods required for the implementation of KeyListener
		public void keyTyped(KeyEvent e) {
		} // NOP
		
		public void keyPressed(KeyEvent e) {
			dispatchKey(e);
		}
		
		public void keyReleased(KeyEvent e) {
		} // NOP
		
		/**
		 * Handles keyboard events, e.g., spacebar toggles
		 * simulation/pausing, and escape resets the current Task.
		 */
		public void dispatchKey(KeyEvent e) {
			switch (e.getKeyCode()) {
				case KeyEvent.VK_S:
					single = true;
					break;
				case KeyEvent.VK_SPACE:
					simulate = !simulate;
					if (simulate) {
						System.out.println("Starting simulation...");
					}
					else {
						System.out.println("Simulation paused.");
					}
					break;
				case KeyEvent.VK_ESCAPE:
					taskSelector.resetToRest(); //sets task=null;
					break;
				case KeyEvent.VK_E:
					frameExporter = ((frameExporter == null) ? (new FrameExporter()) : null);
					System.out.println("'e' : frameExporter = " + frameExporter);
					break;
				case KeyEvent.VK_I:
					frameExporter = ((frameExporter == null) ? (new FrameExporter(true)) : null);
					
					System.out.println("'i' : frameExporter = " + frameExporter);
					break;
				case KeyEvent.VK_L:
					loadFrameFromFile();
					break;
				case KeyEvent.VK_EQUALS:
					N_STEPS_PER_FRAME = Math.max((int) (1.05 * N_STEPS_PER_FRAME), N_STEPS_PER_FRAME + 1);
					System.out.println("N_STEPS_PER_FRAME=" + N_STEPS_PER_FRAME + ";  dt="
						+ (DT / (double) N_STEPS_PER_FRAME));
					break;
				case KeyEvent.VK_MINUS:
					int n = Math.min((int) (0.95 * N_STEPS_PER_FRAME), N_STEPS_PER_FRAME - 1);
					N_STEPS_PER_FRAME = Math.max(1, n);
					System.out.println("N_STEPS_PER_FRAME=" + N_STEPS_PER_FRAME + ";  dt="
						+ (DT / (double) N_STEPS_PER_FRAME));
					break;
				case KeyEvent.VK_LEFT:
					Vector2d vec = new Vector2d(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
					eyePos.x = vec.x * Constants.CAM_COS_THETA - vec.y * Constants.CAM_SIN_THETA + targetPos.x;
					eyePos.z = vec.x * Constants.CAM_SIN_THETA + vec.y * Constants.CAM_COS_THETA + targetPos.z;
					break;
				case KeyEvent.VK_RIGHT:
					vec = new Vector2d(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
					eyePos.x = vec.x * Constants.CAM_COS_THETA + vec.y * Constants.CAM_SIN_THETA + targetPos.x;
					eyePos.z = -vec.x * Constants.CAM_SIN_THETA + vec.y * Constants.CAM_COS_THETA + targetPos.z;
					break;
				
				// TODO(Optional): Make the camera orbit rather than translate?
				case KeyEvent.VK_UP:
					eyePos.y += 1;
					break;
				case KeyEvent.VK_DOWN:
					eyePos.y -= 1;
					break;
				
				default:
			}
		}
		
		/**
		 * "Task" command base-class extended to support
		 * building/interaction via mouse interface.  All objects
		 * extending Task are implemented here as inner classes for
		 * simplicity.
		 *
		 * Add tasks as necessary for different interaction modes.
		 */
		abstract class Task implements MouseListener, MouseMotionListener {
			/** Displays any task-specific OpengGL information,
			 * e.g., highlights, etc. */
			public void display(GL2 gl) {
			}
			
			// Methods required for the implementation of MouseListener
			public void mouseEntered(MouseEvent e) {
			}
			
			public void mouseExited(MouseEvent e) {
			}
			
			public void mousePressed(MouseEvent e) {
			}
			
			public void mouseReleased(MouseEvent e) {
			}
			
			public void mouseClicked(MouseEvent e) {
			}
			
			// Methods required for the implementation of MouseMotionListener
			public void mouseDragged(MouseEvent e) {
			}
			
			public void mouseMoved(MouseEvent e) {
			}
			
			/** Override to specify reset behavior during "escape" button
			 * events, etc. */
			abstract void reset();
			
		}
		
		/** Clicking task that creates particles. */
		class CreateParticleTask extends Task {
			public void mousePressed(MouseEvent e) {
				// TODO(Optional): get the mouse position instead of a random position
				java.util.Random r = new java.util.Random();
				Vector3d x0 = new Vector3d(r.nextFloat(), r.nextFloat(), r.nextFloat());
				Particle lastCreatedParticle = PS.createParticle(x0);
			}
			
			void reset() {
				taskSelector.resetToRest(); //sets task=null;
			}
		}
		
	}
	
	/**
	 * Displays a filechooser, and then loads a frame file.
	 * Files are expected to be in the same format as those exported by the
	 * FrameExporter class.
	 */
	private void loadFrameFromFile() {
		JFileChooser fc = new JFileChooser("./frames");
		int choice = fc.showOpenDialog(frame);
		if (choice != JFileChooser.APPROVE_OPTION) return;
		String fileName = fc.getSelectedFile().getAbsolutePath();
		
		java.io.File file = new java.io.File(fileName);
		if (!file.exists()) {
			System.err.println("Error: Tried to load a frame from a non-existant file.");
			return;
		}
		
		try {
			java.util.Scanner s = new java.util.Scanner(file);
			int numParticles = s.nextInt();
			PS.reset();
			PS.P.clear();
			for (int i = 0; i < numParticles; i++) {
				double x = s.nextDouble();
				double y = s.nextDouble();
				double z = s.nextDouble();
				PS.createParticle(new Vector3d(x, y, z));
			}
			s.close();
			
		}
		catch (Exception e) {
			e.printStackTrace();
			System.err.println("OOPS: " + e);
		}
	}
	
	/**
	 * A class that either writes the current position of all particles to a text file,
	 * or outputs a png of the current window. Toggle the image boolean to switch modes.
	 *
	 * Text file specification:
	 * The file's first line is an integer N denoting the number of particles in the system.
	 * N lines follow, each with 3 floating point numbers describing the points'
	 * x, y, and z coordinates.
	 *
	 * WARNING: the directory "./frames/" must exist for this class to work properly.
	 */
	private class FrameExporter {
		public boolean image = false;
		private int nFrames = 0;
		private long expN = 0;
		private String folder = "";
		
		FrameExporter() {
			expN = System.currentTimeMillis();
			folder = "exports/e" + expN + "/";
			
			File f = new File(folder);
			if (!f.isDirectory()) {
				f.mkdirs();
			}
		}
		
		FrameExporter(boolean image) {
			super();
			this.image = image;
		}
		
		void writeFrame(GL2 gl) {
			long timeNS = System.currentTimeMillis();
			String number = Utils.getPaddedNumber(nFrames, 5, "0");
			String filename = folder + number + (image ? ".png" : ".txt");
			
			try {
				File file = new File(filename);
				if (file.exists()) System.out.println("WARNING: OVERWRITING PREVIOUS FILE: " + filename);
				
				if (image) {
					GLReadBufferUtil rbu = new GLReadBufferUtil(false, false);
					rbu.readPixels(gl, false);
					rbu.write(file);
				}
				else {
					java.io.BufferedWriter output = new java.io.BufferedWriter(new java.io.FileWriter(file));
					
					output.write("" + PS.P.size() + "\n");
					for (Particle p : PS.P) {
						output.write("" + p.x.x + " " + p.x.y + " " + p.x.z + "\n");
					}
					output.close();
				}
				
				System.out.println((timeNS - expN) + "ms:  Wrote frame: " + filename);
				
			}
			catch (Exception e) {
				e.printStackTrace();
				System.out.println("OOPS: " + e);
			}
			
			nFrames += 1;
		}
	}
	
	/**
	 * ### Runs the ParticleSystemBuilder. ###
	 */
	public static void main(String[] args) {
		try {
			ParticleSystemBuilder psb = new ParticleSystemBuilder();
			psb.start();
			
		}
		catch (Exception e) {
			e.printStackTrace();
			System.out.println("OOPS: " + e);
		}
	}
}
