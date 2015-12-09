#!/bin/sh

for i in $(ls particles_*); do
    python visualizer.py $i
done

ffmpeg -framerate 30 -i frame_%05d.png -vcodec libx264 -r 30 -pix_fmt yuv420p -y movie.mp4

