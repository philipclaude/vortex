#!/bin/sh

vortex=$1

# create adapted triangle mesh of oceans and continents
if [ ! -f data/earth.meshb ]; then
  # pass 0.02 in workflow, but 0.005 for better resolution
  $vortex mesh data/oceans_2048.png --hmin $2 --output data/earth.meshb
fi

# extract oceans and land
$vortex extract data/earth.meshb --oceans data/oceans.meshb --continents data/land.meshb

# voronoi diagrams
$vortex voronoi --domain sphere --points random --n_points 100000 --n_smooth 10 --output data/voronoi_sphere.meshb
$vortex voronoi --domain data/oceans.meshb --points random --n_points 100000 --n_smooth 10 --output data/voronoi_oceans.meshb
$vortex voronoi --domain data/oceans.meshb --points random --n_points 1000 --n_smooth 10 --output data/voronoi_oceans_coarse.meshb

# merge
$vortex merge data/voronoi_oceans.meshb --combine --output data/voronoi_oceans_merged.meshb
$vortex merge data/voronoi_oceans_coarse.meshb --combine --output data/voronoi_oceans_coarse_merged.meshb

