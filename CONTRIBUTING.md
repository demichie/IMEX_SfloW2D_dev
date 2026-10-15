# Contributing to IMEX_SfloW2D

Thank you for your interest in improving IMEX_SfloW2D, a depth-averaged numerical model for geophysical and volcanic flows. Contributions to the Fortran solver, tests, reproducible examples, scientific documentation, and development tools are welcome.

## Before contributing

- Search the [existing issues](https://github.com/demichie/IMEX_SfloW2D_dev/issues) and pull requests for related work.
- For significant numerical, physical, or architectural changes, open an issue first to discuss the mathematical formulation, expected behavior, validation strategy, and compatibility implications.
- Keep changes focused. Separate file relocations, refactoring without intended numerical changes, new numerical methods, and scientific model extensions into distinct pull requests where practical.
- Never commit credentials, machine-specific paths, generated simulation outputs, or large data files unless they are intentional, documented, and required fixtures.

## Report a bug

[Open an issue](https://github.com/demichie/IMEX_SfloW2D_dev/issues/new) and include the following where applicable:

- The source revision or release, operating system, compiler and version, build profile, and relevant library versions (especially NetCDF-Fortran).
- A minimal input file or existing example that reproduces the problem, with the steps required to prepare it.
- The exact command, expected behavior, observed behavior, and complete relevant error output.
- The OpenMP thread count and whether the problem occurs in a restart or in a continuous run.
- For numerical problems, the affected physical quantities, time interval, and an appropriate reference result or conservation diagnostic.

Do not attach confidential information or prohibitively large simulation outputs; describe how maintainers can reproduce the issue instead.

## Build and test

The project uses Autotools and requires a compatible Fortran compiler, NetCDF-Fortran, LAPACK/BLAS, and the usual build utilities. From the repository root, a representative GNU Fortran build on Linux is:

```sh
autoreconf -fi
./configure --with-netcdf=/usr FC=gfortran
make -j2
make check
```

For a diagnostic GNU build, use a clean build directory or clean the preceding build and configure with:

```sh
./configure --with-netcdf=/usr FC=gfortran --enable-strict-debug
```

`FCFLAGS`, `OPENMP_FCFLAGS`, `NETCDF_LIBS`, and `LAPACK_LIBS` can be supplied for a different toolchain. Consult [README.md](README.md) and the CI workflows for additional details. A NetCDF-Fortran module built with one Fortran compiler may not be compatible with another.

Relevant executable regression suites include:

```sh
TESTS/run_tests.sh src/IMEX_SfloW2D
TESTS/TEST_BINARY_RESTART/run_test.sh src/IMEX_SfloW2D
TESTS/TEST_STOCHASTIC_RESTART/run_test.sh src/IMEX_SfloW2D
TESTS/TEST_PCCU_GATE_H/run_test.sh src/IMEX_SfloW2D
```

These are in addition to `make check`; they do not replace tests targeted at the code being modified. The GitHub Actions workflows validate GNU optimized/strict-debug builds, Docker configurations, and Intel `ifx` portability within their recorded scopes. These validation workflows run on tag pushes, published releases, or manual dispatch; they do not run on every branch commit or pull request.

## Numerical changes and scientific validation

IMEX_SfloW2D maintains a documented scientific baseline and acceptance evidence under [`docs/acceptance/`](docs/acceptance/). Contributions affecting the solver must:

1. Explain whether a change is intended to preserve numerical results or deliberately alter the equations, discretization, source terms, boundary conditions, or solver algorithm.
2. Include focused regression tests covering the changed functionality, including relevant one- and two-dimensional cases, dry/wet interfaces, conservation, restart continuity, and OpenMP behavior when applicable.
3. Preserve existing supported input behavior by default, or explicitly document and test any intentional compatibility change.
4. Compare relevant results to the accepted reference revision with the appropriate scientific tolerances; do not assume bitwise identity across compilers or optimization modes.
5. Record test conditions, reference revisions, deviations, and any new limitations. A passing build alone does not constitute scientific acceptance.

Pure source relocation should avoid editing Fortran formulas or changing execution order in the same commit. Performance-oriented refactoring should also assess extra allocations, array temporaries, and runtime cost.

## Fortran source documentation

Document new or changed public modules, types, procedures, and interfaces using the repository's Doxygen conventions. Include the purpose, physical or mathematical meaning, units where relevant, and the role and intent of procedure arguments (including type-bound `this` arguments). Update documentation together with the code rather than deferring it to a later change.

The source reference is generated from [`template.doxy`](template.doxy) and published with the documentation website. The documentation build and the source-documentation tests complement, but do not replace, scientific validation.

## Adding or changing an example

**Every new runnable example must include an `EXAMPLES/<example-path>/README.md`** following the repository's standardized Markdown format. This applies to nested cases as well as top-level example directories. The documentation website scans example directories recursively and generates the public Examples catalogue directly from these files.

Start the file with YAML front matter containing the required fields:

```yaml
---
title: "Descriptive example title"
summary: "One- or two-sentence summary of the simulation."
category: "Numerical benchmarks"
tags:
  - "validation"
  - "shallow flow"
---
```

Then include these sections, with the exact headings shown below so that documentation tools can recognize them:

- `## Overview` — physical problem, objective, and relevant assumptions.
- `## Main parameters` — a readable table of the key physical and numerical settings, with units and values consistent with the actual input or generator scripts.
- `## How to run` — required files, dependencies, preprocessing, commands to run from a clearly stated directory, and any necessary postprocessing.
- `## Expected output` — generated files, their formats, and the main results or behavior to examine.

Add `## References` whenever the case derives from published work or cites a scientific source. Retain additional sections, such as physical and numerical details, validation, limitations, and postprocessing, whenever they contain useful information. **Do not discard existing bibliographic references, experimental context, parameter explanations, or operational instructions when updating a README.** Preserve existing `README.txt` documentation until its complete content has been checked against the standardized replacement.

Describe what the supplied files actually do; do not infer undocumented parameter values from directory names or present unverified numerical results as validated. Where scripts generate input or DEM data, explain preparation and dependencies sufficiently for another user to reproduce the example.

For an existing example, update its `README.md` whenever changing its inputs, preprocessing, scripts, outputs, or scientifically significant defaults. Consult [`EXAMPLES/README_SECOND_PASS_AUDIT.md`](EXAMPLES/README_SECOND_PASS_AUDIT.md) for outstanding historical documentation checks.

## Documentation website

The website uses MkDocs Material and integrates Doxygen. The Examples pages are generated from `EXAMPLES/**/README.md` rather than maintained as duplicate copies. Static website content is located under `docs/website/`, and the generator is under `tools/docs/`.

Publication is performed by [the documentation workflow](.github/workflows/publish-documentation.yml), which can be run manually or on tag/release events; it is **not** run for every ordinary commit. When modifying the generator, metadata rules, or page layout, verify the generated Examples catalogue and documentation build. Do not publish output directories or generated HTML directly into source folders.

## Pull requests

Use a descriptive title and explain:

- The problem or feature addressed, and the affected modules or examples.
- Whether numerical output is expected to change, and why.
- The tests performed, toolchain information, and relevant evidence or result comparisons.
- Any modifications to input formats, dependencies, public documentation, or restart compatibility.
- Any known limitations or follow-up work.

Keep unrelated formatting changes out of numerical patches. Wait for applicable CI checks and address failures before requesting merge. Maintainers may require additional acceptance tests for changes that affect numerical or scientific behavior.

## License and attribution

The repository contains licensing and authorship information in [`COPYING`](COPYING), [`AUTHORS`](AUTHORS), and [`CITATION.cff`](CITATION.cff). Preserve copyright and attribution notices and cite scientific references relevant to contributed algorithms and examples.
