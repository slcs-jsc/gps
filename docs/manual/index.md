# GPS Code Collection

The GPS Code Collection provides software for the analysis of GPS radio occultation observations.

This MkDocs manual is the lightweight user-facing entry point for the repository. It complements:

- the top-level `README.md` for installation notes
- the generated Doxygen reference for implementation details
- the regression tests in `tests/` for concrete workflow examples

## What You Will Find Here

- a quick build and test guide
- an overview of the available command-line tools
- links to the main project resources

## Main Components

The repository is organized around a C-based command-line toolkit in `src/`, together with:

- local dependency builds in `libs/`
- regression tests in `tests/`
- documentation sources in `docs/`

The source tree builds several executables for working with profile data, maps, perturbations, tropopause diagnostics, and related statistics.

## External Components

The files `src/jurassic.c` and `src/jurassic.h` are included here for convenience from the external [JURASSIC radiative transfer model](https://github.com/slcs-jsc/jurassic).

## Contact

We are interested in sharing the GPS Code Collection for research applications. For questions or support, please contact:

Dr. Lars Hoffmann  
Jülich Supercomputing Centre, Forschungszentrum Jülich  
<l.hoffmann@fz-juelich.de>
