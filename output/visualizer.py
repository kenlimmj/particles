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

def main(fname):
    if fname[0:15] != "data/particles_":
        return
    frame = fname[15:20]
    outname = "frame/" + frame + ".png"

    fin = open(fname, "r")
    n = int(fin.readline().strip())
    bs = float(fin.readline().strip())
    px = np.ndarray(n)
    py = np.ndarray(n)
    pz = np.ndarray(n)
    count = 0
    for l in fin.readlines():
        parts = l.strip().split(" ")
        px[count] = float(parts[0])
        py[count] = float(parts[1])
        pz[count] = float(parts[2])
        count = count + 1
    fin.close()

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-0.25 * bs, 1.25 * bs)
    ax.set_ylim(-0.25 * bs, 1.25 * bs)
    ax.set_zlim(0.0, bs * 1.5)
    ax.scatter(px, pz, py, s=10*10, c=[(0, 0, 1, 0.5)] * n)
    fig.savefig(outname)


if len(sys.argv) > 1:
    main(sys.argv[1])

