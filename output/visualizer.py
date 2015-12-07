#!/usr/bin/env python

"""
Visualize SPH.

NB: Requires a modern Matplotlib version; also needs
 either FFMPeg (for MP4) or ImageMagick (for GIF)
"""
import glob, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation
import sys
    
# For each data file in ./Data, output heights using data2height (an NGA binary)
times = list()
f = 0
for files in os.listdir("."):
  if files[0:10] == "particles_":
    times.append(files[10:15])
    f = f + 1

with open("./particles_"+times[0],'r') as p:
  n = p.readline().strip()
n = int(n)
locations = np.ndarray((n,3))
velocities = np.ndarray((n,3))

fig = plt.figure(figsize=(10,10))
def plot_frame(i):
  ax = fig.add_subplot(111, projection='3d')
  ax.set_xlim(0.0, 1.0)
  ax.set_ylim(0.0, 1.0)
  ax.set_zlim(0.0, 2.0)
  data = np.loadtxt("./particles_"+times[i], skiprows=1)
  locations[:,0] = data[:,0]
  locations[:,1] = data[:,1]
  locations[:,2] = data[:,2]
  #velocities[:,0] = data[:,3]
  #velocities[:,1] = data[:,4]
  #velocities[:,2] = data[:,5]
  X = data[:,0]
  Y = data[:,1]
  Z = data[:,2]
  ax.scatter(X, Z, Y, s=10*10)
  return ax

metadata = dict(title='SPH', artist='Matplotlib')

Writer = manimation.writers['ffmpeg']
writer = Writer(fps=15, metadata=metadata,
                extra_args=["-r", "30",
                            "-c:v", "libx264",
                            "-pix_fmt", "yuv420p"])

with writer.saving(fig, "movie.mp4", f):
  for i in range(f):
    print "frame ",i," written. ",float(i)/float(f)*100.0,"% done. \n",
    ax = plot_frame(i)
    writer.grab_frame()
    plt.delaxes(ax)
