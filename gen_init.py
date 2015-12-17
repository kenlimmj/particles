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
nprocx = 2               # processors in x
nprocy = 2               # processors in y
nprocz = 2               # processors in z
boxsize = 1.0
nproc=nprocx*nprocy*nprocz

######## NOTHING NEEDS TO BE CHANGED AFTER THIS POINT #######
pdist = h/1.3
nx = int(math.ceil(npart**(1.0/3.0)))
ny = nx
nz = nx
n = nx*ny*nz   # Total number of particles

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

px = np.linspace(Lix,Lfx,nx)
py = np.linspace(Liy,Lfy,ny)
pz = np.linspace(Liz,Lfz,nz)
vx[:] = Ivx
vy[:] = Ivy
vz[:] = Ivz

n_local=[0] * (nproc)

for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
        proc = (nprocx*nprocy*int((pz[k]/(boxsize/nprocz)))
             +  nprocx*       int((py[j]/(boxsize/nprocy)))
             +                int((px[i]/(boxsize/nprocx))))
        n_local[proc] = n_local[proc]+1

# Write out to file
fid=[]
for proc in range(0,nproc):
  procn='00000'
  temp=str(proc)
  procn=procn[0:len(procn)-len(temp)]+temp
  fid.append(open(filename+procn, 'w'))
  fid[proc].write(str(n))
  fid[proc].write("\n")
  fid[proc].write(str(n_local[proc]))
  fid[proc].write("\n")
  fid[proc].write(str(h))
  fid[proc].write("\n")
  fid[proc].write(str(tfin))
  fid[proc].write("\n")
  fid[proc].write(str(dt))
  fid[proc].write("\n")
  fid[proc].write(str(ff))
  fid[proc].write("\n")
  fid[proc].write(str(nprocx))
  fid[proc].write("\n")
  fid[proc].write(str(nprocy))
  fid[proc].write("\n")
  fid[proc].write(str(nprocz))
  fid[proc].write("\n")
  fid[proc].write(str(boxsize))
  fid[proc].write("\n")
for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
        proc = (nprocx*nprocy*int((pz[k]/(boxsize/nprocz)))
             +  nprocx*       int((py[j]/(boxsize/nprocy)))
             +                int((px[i]/(boxsize/nprocx))))
        line = str(px[i])+" "+str(py[j])+" "+str(pz[k])+" "+str(vx[i])+" "+str(vy[j])+" "+str(vz[k])
        fid[proc].write(line)
        fid[proc].write("\n")
for proc in range(0,nproc):
  fid[proc].close()

