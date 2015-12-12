#!/usr/bin/env python

# Configure MPL. This has to appear before anything else
# to ensure that MPL does not try to use Xwindows to render
import matplotlib
matplotlib.use('Agg')

import matplotlib.animation as manimation
import matplotlib.pyplot as plt

from datetime import datetime
from glob import glob
from mpl_toolkits.mplot3d import Axes3D
from multiprocessing import Pool
from numpy import ndarray
from os import getcwd
from os.path import join
from subprocess import call

# The number of iteration steps taken
# NUM_ITER_STEPS = 100

# Directory paths
OUTPUT_DIR = join(getcwd(), 'output')
DATA_DIR = join(OUTPUT_DIR, 'data')

# File names
OUTPUT_FNAME = 'frame_{}.png'
DATA_FNAME_GLOB = 'particles_*'

def procParticleData(f):
    prefix, suffix = f.split('_')
    outFileName = OUTPUT_FNAME.format(suffix)

    with open(f, 'r') as fid:
        n = int(fid.readline().strip())
        px, py, pz = ndarray(n), ndarray(n), ndarray(n)
        count = 0

        for l in fid.readlines():
            parts = l.strip().split()
            px[count], py[count], pz[count] = float(parts[0]), float(parts[1]), float(parts[2])
            count += 1

    fig = plt.figure(figsize=(10, 10))
    ax= fig.add_subplot(111, projection='3d')

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_zlim(0.0, 2.0)

    ax.scatter(px, py, pz, s=10*10, c=[(0, 0, 1, 0.5)]*n)
    fig.savefig(join(DATA_DIR, outFileName))
    plt.close(fig)


def main():
    # Extract all the particle data files
    particleFiles = glob(join(DATA_DIR, DATA_FNAME_GLOB))

    # The number of particle data files should be the number of iterations performed
    # assert(len(particleFiles) == NUM_ITER_STEPS + 1)

    # Generate figures in parallel BECAUSE WE BE PARALLEL PROGRAMMIN' YO'
    Pool().map(procParticleData, particleFiles)

    dt = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')

    # Generate a movie from the still frames
    call([
        'ffmpeg',
        '-framerate', '100',
        '-i', join(DATA_DIR, 'frame_%05d.png'),
        '-vcodec', 'libx264',
        '-r', '60',
        '-pix_fmt', 'yuv420p',
        '-y', join(OUTPUT_DIR, 'movie_{}.mp4'.format(dt))
    ])


if __name__ == '__main__':
    main()
