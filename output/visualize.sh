#!/bin/sh

#for i in $(ls data/particles_*); do
#    python visualizer.py $i
#done
ls data/particles_* | xargs -n 1 -P 24 python visualizer.py

ffmpeg -framerate 30 -i frame/%05d.png -vcodec libx264 -r 30 -pix_fmt yuv420p -y movie.mp4

