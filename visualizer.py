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
    times.append(files[10:19])
    f = f + 1

tmp = np.loadtxt("./particles_"+times[0],skiprows=1)
n = len(tmp[:,0])
locations = np.ndarray((n,3,f))
velocities = np.ndarray((n,3,f))
for t in range(0,f):
  data = np.loadtxt("./particles_"+times[t], skiprows=1)
  locations[:,0,t] = data[:,0]
  locations[:,1,t] = data[:,1]
  locations[:,2,t] = data[:,2]
  velocities[:,0,t] = data[:,3]
  velocities[:,1,t] = data[:,4]
  velocities[:,2,t] = data[:,5]


fig = plt.figure(figsize=(10,10))

def plot_frame(i, stride=1):
  ax = fig.add_subplot(111, projection='3d')
  ax.set_zlim(0, 0.1)
  X = locations[:,0,i]
  Y = locations[:,1,i]
  Z = locations[:,2,i]
  ax.scatter(X, Z, Y)
  return ax

metadata = dict(title='SPH', artist='Matplotlib')

Writer = manimation.writers['ffmpeg']
writer = Writer(fps=15, metadata=metadata,
                extra_args=["-r", "30",
                            "-c:v", "libx264",
                            "-pix_fmt", "yuv420p"])

with writer.saving(fig, "movie.mp4", f):
  for i in range(f):
    ax = plot_frame(i)
    writer.grab_frame()
    plt.delaxes(ax)
