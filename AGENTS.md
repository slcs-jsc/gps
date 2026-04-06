# AGENTS.md

Guidance for coding agents working in this repository.

## Scope

This file applies to the entire repository rooted here.

## Repository Overview

This project is the GPS Code Collection, a C codebase for analysis of GPS/RO observations.

Important directories:

- `src/`: primary C sources, headers, and the main `Makefile`
- `tests/`: shell-based regression tests and reference outputs
- `docs/`: Doxygen and MkDocs documentation sources
- `libs/`: vendored third-party source archives and the local dependency build script

Important external source files:

- `src/jurassic.c`
- `src/jurassic.h`

These two files are copied into this repository for convenience from the external JURASSIC radiative transfer model hosted on GitHub. Treat them as vendored upstream code. Avoid editing them unless the task explicitly requires syncing or patching that external dependency, and keep any such changes clearly isolated and documented.
Upstream repository: `https://github.com/slcs-jsc/jurassic`

## Build Workflow

Primary build commands are run from `src/`:

```sh
cd src
make
make check
```

Useful build flags supported by `src/Makefile`:

- `STATIC=1`: enable static linking
- `OPT=...`: override optimization level
- `INFO=1`: emit optimization info
- `PROF=1`: enable profiling flags
- `COV=1`: enable coverage instrumentation

The executables are built in `src/` and can be installed to `../bin` via:

```sh
cd src
make install
```

## Tests

Regression tests are orchestrated by `src/Makefile` and executed from `tests/`:

- `make check` runs `pert_test` and `tropo_test`
- Each test has a `run.sh` script under `tests/<name>/`

Test scripts currently:

- set `LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH`
- set `OMP_NUM_THREADS=4`
- create a local `data/` directory inside the test folder
- compare generated outputs against `data.ref/`

If you change numerical behavior or output formatting, run the relevant regression tests and check diffs carefully before updating any reference data.

## Dependencies

This project expects GSL and netCDF headers/libraries. The repository includes vendored archives in `libs/` and a helper script:

```sh
cd libs
./build.sh
```

By default, `src/Makefile` looks for headers and libraries in:

- `../libs/build/include`
- `../libs/build/lib`
- `../libs/build/lib64`

Do not modify vendored archives or rebuilt library artifacts unless the task is explicitly about dependency updates or build system changes.
Apply the same caution to vendored upstream source files such as `src/jurassic.c` and `src/jurassic.h`.

## Code Style

Follow the existing style in `src/`:

- C, not C++
- keep functions and declarations in the project’s current formatting style
- use the shared helpers/macros already present in headers such as `NC`, `ALLOC`, and error-reporting utilities where appropriate
- prefer minimal, localized changes over large refactors
- preserve Doxygen-style comments and block comment conventions when touching documented code

The `Makefile` uses strict warnings and `-Werror`. New code should be clean under those flags, including conversion and prototype warnings.

## File Handling

- Treat `tests/data/` and `tests/*/data.ref/` as regression fixtures
- Avoid unnecessary edits to generated binaries, archives, PDFs, or large netCDF reference files
- Keep changes focused on source, build, test, or docs files unless the task explicitly requires fixture updates

## Documentation

Documentation tooling is driven from `src/Makefile`:

- `make doxygen`
- `make mkdocs`

If code changes affect CLI behavior, build instructions, or documented algorithms, update the relevant docs in `README.md` or `docs/`.

## Agent Expectations

- Start by reading `README.md` and `src/Makefile` for build context
- Prefer `rg` for code search
- Before changing behavior, identify which executable in `src/` owns it
- Validate with the smallest relevant command first, then broaden to `make check` when feasible
- Do not revert unrelated user changes in the worktree
