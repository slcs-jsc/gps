# Tools

The `src/Makefile` builds the following command-line executables:

- `day2doy`
- `doy2day`
- `events`
- `jsec2time`
- `map`
- `perturbation`
- `prof`
- `response`
- `time2jsec`
- `tropo`
- `variance`

## What They Are Used For

The tools cover a small set of common workflows around GPS radio occultation data:

- date and time conversion utilities
- event extraction and tabular summaries
- map and profile export
- perturbation analysis
- tropopause diagnostics
- response and variance calculations

## Related Components

Several executables share infrastructure from:

- `libgps.c` and `libgps.h`: project-specific I/O and analysis helpers
- `jurassic.c` and `jurassic.h`: vendored code from the external [JURASSIC](https://github.com/slcs-jsc/jurassic) radiative transfer model

## Build And Test Entry Points

Most users will interact with the tools through the main build workflow:

```sh
cd src
make
make check
```

For implementation details and APIs, see the Doxygen manual linked from the documentation site.
