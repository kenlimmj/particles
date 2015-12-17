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
filename = sys.argv[1]   # File name to write
box = 5.0                # Box size
Lix = 0.0                # Starting of fluid in x
Liy = 0.0                # Starting of fluid in y
Liz = 0.0                # Starting of fluid in z
Lfx = 2.18               # End of fluid in x
Lfy = 2.18               # End of fluid in y
Lfz = 2.18               # End of fluid in z
Ivx = 0.0                # Initial x velocity
Ivy = 0.0                # Initial y Velocity
Ivz = 0.0                # Initial z velocity



######## NOTHING NEEDS TO BE CHANGED AFTER THIS POINT #######
pdist = 0.06
Lfx = Lfx - pdist
Lfy = Lfy - pdist
Lfz = Lfz - pdist
nx = int((Lfx-Lix)/pdist + 0.5) + 1
ny = int((Lfy-Liy)/pdist + 0.5) + 1
nz = int((Lfz-Liz)/pdist + 0.5) + 1
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
fid.write(str(box))
fid.write("\n")
for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
      line = str(px[i])+ " "+str(py[j])+" "+str(pz[k])+" "+str(vx[i])+" "+str(vy[j])+" "+str(vz[k])
      fid.write(line)
      fid.write("\n")
fid.close()

