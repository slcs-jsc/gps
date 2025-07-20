#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
gps=../../src

# Create directory...
rm -rf data && mkdir -p data

# Create perturbation file...
$gps/tropo - data/tropo.nc ../data/*_nc

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.nc) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
