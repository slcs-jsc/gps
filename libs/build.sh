#! /bin/bash

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.gz $target/src \
    && cd $target/src \
    && for f in $(ls *tar.gz) ; do tar xvzf $f ; done \
	|| exit

# GSL...
dir=gsl-2.7.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# netCDF...
dir=netcdf-c-4.6.2
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-dap --disable-netcdf-4 \
    && make -j && make check && make install && make clean \
	|| exit
