at#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
gps=../../src

# Create directory...
rm -rf data && mkdir -p data

# Create perturbation file...
$gps/perturbation - data/pert.nc ../data/*_nc METBASE - HAM_DZ 10 HAM_DZ2 1

# Extract map...
$gps/map - data.ref/pert.nc data/map.tab

# Extract profile...
$gps/prof - data/pert.nc data/prof.tab IDS 10

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.nc data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
