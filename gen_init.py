#!/usr/bin/env python

import sys
import string
import numpy as np
import shutil
import glob, os
import subprocess
import matplotlib.pyplot as plt
import math

h = 0.05
pdist = h/1.3
Lx = 0.5   # Particles along x
Ly = 0.5  # Particles along y
Lz = 1  # Particles along z
nx = int(math.ceil(float(Lx)/pdist))
ny = int(math.ceil(float(Ly)/pdist))
nz = int(math.ceil(float(Lz)/pdist))

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

px = np.linspace(0.0+pdist,Lx,nx)
py = np.linspace(0.0+pdist,Ly,ny)
pz = np.linspace(0.0+pdist,Lz,nz)
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

