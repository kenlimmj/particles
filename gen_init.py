import sys
import string
import numpy as np
import shutil
import glob, os
import subprocess
import matplotlib.pyplot as plt

nx = 25   # Particles along x
ny = 25   # Particles along y
nz = 25  # Particles along z
n = nx*ny*nz   # Total number of particles

BS = 1.0   # Box Size
buffer = 0.125

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

px = np.linspace(-BS+buffer,BS-buffer,nx)
py = np.linspace(0+buffer,2*BS-buffer,ny)
pz = np.linspace(-BS+buffer,BS-buffer,nz)
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

