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

Run `make clean` before rebuilding with changed compiler options; Make does not track command-line flag changes. Run `make check` without `-j`.

## Tests

Regression tests are orchestrated by `src/Makefile` and executed from `tests/`:

- `make check` runs `pert_test` and `tropo_test`
- Each test has a `run.sh` script under `tests/<name>/`

Test scripts currently:

- set `LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH`
- set `OMP_NUM_THREADS=4`
- create a local `data/` directory inside the test folder
- compare generated outputs against `data.ref/`

If you change numerical behavior or output formatting, run the relevant regression tests and investigate diffs before updating any reference data. Update only affected baselines after verifying and explaining the new results.

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
- before changing behavior, check existing callers, tests, and documented expectations; question assumptions and explain non-obvious choices
- reuse existing data structures where possible; check memory and affected data paths before adding a `struct`
- avoid new `static` helper functions; use a macro for a small operation when its arguments can be evaluated safely, or a regular function declared and documented in `src/libgps.h`
- preserve Doxygen-style comments and block comment conventions when touching documented code
- format C sources and headers from `src/` with `make indent` after a successful build; review the diff because the target rewrites all `*.c` and `*.h` files, including copied JURASSIC sources

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
Keep the MkDocs manual in `docs/manual/` and Doxygen API documentation current with relevant code, build, and workflow changes. Build affected docs with `make mkdocs` and/or `make doxygen` from `src/` when the tools are available. For new or revised API Doxygen comments, use a substantive description beyond `@brief`, document parameters and behavior where relevant, and include `@author` with the responsible author's name. Follow the detailed `/** ... */` examples in `src/jurassic.h` without editing that copied dependency for style alone; do not invent attribution.

## Agent Expectations

- Start by reading `README.md` and `src/Makefile` for build context
- Prefer `rg` for code search
- Before changing behavior, identify which executable in `src/` owns it
- Validate with the smallest relevant command first, then broaden to `make check` when feasible
- Do not revert unrelated user changes in the worktree
