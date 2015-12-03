#!/usr/bin/env python

import sys
import string
import numpy as np
import shutil
import glob, os
import subprocess
import matplotlib.pyplot as plt

nx = 300   # Particles along x
ny = 1  # Particles along y
nz = 1  # Particles along z
n = nx*ny*nz   # Total number of particles

BS = 0.5   # Box Size
buffer = 1.0e-5

# File name to write
filename = "init.txt"

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

px = np.linspace(-BS+buffer,0.0,nx)
py[:] = 0.0
#py = np.linspace(0+buffer,BS-buffer,ny)
#pz = np.linspace(-BS+buffer,BS-buffer,nz)
pz[:] = 0.0
vx[:] = 0.0
vy[:] = 0.0
vz[:] = 0.0

# Write out to file
fid = open(filename, 'w')
fid.write(str(n))
fid.write("\n")
for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
      line = str(px[i])+ " "+str(py[j])+" "+str(pz[k])+" "+str(vx[i])+" "+str(vy[j])+" "+str(vz[k])
      fid.write(line)
      fid.write("\n")
fid.close()

