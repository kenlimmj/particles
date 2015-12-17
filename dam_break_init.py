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
filename = "init"    # File name to write
tfin = 1.00              # End time of simulation
dt = 0.002              # Time step
ff = 0.01                # Time between viz frames
h = 0.05                 # Kernel (particle) size
npart = 35152               # Number of particles


######## NOTHING NEEDS TO BE CHANGED AFTER THIS POINT #######
pdist = h/1.3
nx = int(math.ceil(npart**(1.0/3.0)))
ny = nx
nz = nx
n = nx*ny*nz   # Actual Total number of particles

# Initialize
px = np.ndarray((n))
py = np.ndarray((n))
pz = np.ndarray((n))
vx = np.ndarray((n))
vy = np.ndarray((n))
vz = np.ndarray((n))

Lix = 0.0
Liy = 0.0
Liz = 0.0
Ivx = 0.0
Ivy = 0.0
Ivz = 0.0
Lfx = pdist*float(nx)
Lfy = pdist*float(ny)
Lfz = pdist*float(nz)
bs = Lfx*2

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
fid.write(str(bs))
fid.write("\n")
for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
      line = str(px[i])+ " "+str(py[j])+" "+str(pz[k])+" "+str(vx[i])+" "+str(vy[j])+" "+str(vz[k])
      fid.write(line)
      fid.write("\n")
fid.close()

