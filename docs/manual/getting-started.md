# Getting Started

This page summarizes the usual setup, build, and test workflow for the GPS Code Collection on Linux.

## Repository Layout

The repository is organized into a few main areas:

- `src/`: C source files, headers, and the main `Makefile`
- `tests/`: regression tests and reference outputs
- `libs/`: vendored dependency archives and the local build helper
- `docs/`: Doxygen and MkDocs documentation sources

## Prerequisites

You need a standard C toolchain, including:

- `gcc`
- `make`

The code also depends on:

- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
- [netCDF C library](https://www.unidata.ucar.edu/software/netcdf/)

If those libraries are not already available on your system, the repository provides local copies and a helper build script.

## Build Dependencies

To build the vendored libraries locally:

```sh
cd libs
./build.sh
```

The default source build looks for headers and libraries in:

- `../libs/build/include`
- `../libs/build/lib`
- `../libs/build/lib64`

## Build The Code

Compile the executables from `src/`:

```sh
cd src
make
```

Useful build options from the `Makefile`:

- `STATIC=1`: enable static linking
- `OPT=-O0` or `OPT=-O3`: change optimization level
- `INFO=1`: emit optimization information
- `PROF=1`: build with profiling flags
- `COV=1`: build with coverage instrumentation

Example:

```sh
cd src
make OPT=-O0 COV=1
```

## Run Regression Tests

The repository includes shell-based regression tests:

```sh
cd src
make check
```

This runs the current test suite:

- `pert_test`
- `tropo_test`

You can also run a single test directly:

```sh
cd tests/pert_test
./run.sh
```

The test scripts create temporary output in a local `data/` directory and compare it with checked-in reference results under `data.ref/`.

## Install Executables

The default install target copies executables to `../bin`:

```sh
cd src
make install
```

To override the target directory:

```sh
cd src
make DESTDIR=/path/to/bin install
```

## Notes

- The build uses strict compiler warnings and treats warnings as errors via `-Werror`.
- Static linking may require adjustments on some systems; if it causes issues, review the `STATIC` setting and linker flags in `src/Makefile`.
- The files `src/jurassic.c` and `src/jurassic.h` are vendored from the external [JURASSIC](https://github.com/slcs-jsc/jurassic) radiative transfer model.
