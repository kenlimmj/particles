#!/usr/bin/env python

import sys
import string
import numpy as np
import shutil
import glob, os
import subprocess
import matplotlib.pyplot as plt
import math

######## INPUT PARAMETERS ###############
# File name to write
filename = "init.txt"    # File name to write
tfin = 3.00              # End time of simulation
dt = 0.0001              # Time step
ff = 0.05                # Time between viz frames
h = 0.05                 # Kernel (particle) size
Lix = 0.00               # Starting of fluid in x
Liy = 0.00               # Starting of fluid in y
Liz = 0.00               # Starting of fluid in z
Lfx = 0.50               # End of fluid in x
Lfy = 0.75               # End of fluid in y
Lfz = 0.25               # End of fluid in z
Ivx = 0.0                # Initial x velocity
Ivy = 0.0                # Initial y Velocity
Ivz = 0.0                # Initial z velocity



######## NOTHING NEEDS TO BE CHANGED AFTER THIS POINT #######
pdist = h/1.3
nx = int(math.ceil(float(Lfx-Lix)/pdist))
ny = int(math.ceil(float(Lfy-Liy)/pdist))
nz = int(math.ceil(float(Lfz-Liz)/pdist))
n = nx*ny*nz   # Total number of particles

# Initialize
px = np.ndarray((n))
py = np.ndarray((n))
pz = np.ndarray((n))
vx = np.ndarray((n))
vy = np.ndarray((n))
vz = np.ndarray((n))


##### Dam Break #######
# Start all particles on left half
# Evenly distribute in height and depth

px = np.linspace(Lix,Lfx,nx)
py = np.linspace(Liy,Lfy,ny)
pz = np.linspace(Liz,Lfz,nz)
vx[:] = Ivx
vy[:] = Ivy
vz[:] = Ivz

# Write out to file
fid = open(filename, 'w')
fid.write(str(n))
fid.write("\n")
fid.write(str(h))
fid.write("\n")
fid.write(str(tfin))
fid.write("\n")
fid.write(str(dt))
fid.write("\n")
fid.write(str(ff))
fid.write("\n")
for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
      line = str(px[i])+ " "+str(py[j])+" "+str(pz[k])+" "+str(vx[i])+" "+str(vy[j])+" "+str(vz[k])
      fid.write(line)
      fid.write("\n")
fid.close()

