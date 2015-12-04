#!/usr/bin/env python

import sys
import string
import numpy as np
import shutil
import glob, os
import subprocess
import matplotlib.pyplot as plt
import math

h = 1.0
pdist = h/1.3
Lix = 5
Liy = 5
Liz = 5
Lfx = 7.5
Lfy = 7.5
Lfz = 7.5
nx = int(math.ceil(float(Lfx-Lix)/pdist))
ny = int(math.ceil(float(Lfy-Liy)/pdist))
nz = int(math.ceil(float(Lfz-Liz)/pdist))

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

px = np.linspace(Lix,Lfx,nx)
py = np.linspace(Liy,Lfy,ny)
pz = np.linspace(Liz,Lfz,nz)
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

