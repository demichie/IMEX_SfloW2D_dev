# Manual GitHub Actions

Every workflow in this repository is launched **manually** from the
GitHub **Actions** tab using **Run workflow** (`workflow_dispatch`).
No workflow is triggered by a push, pull request, release, tag or schedule.

## Active workflows

| Workflow | Purpose | When to run |
| --- | --- | --- |
| [Build (optimized and debug)](workflows/build-only.yml) | Compile the full Fortran executable with gfortran in the **selected** `optimized` or `debug` profile (debug enables `--enable-strict-debug`). Exactly one build job runs; **no tests run**. | After source changes, or for a fast compiler check. |
| [Solver integration and OpenMP validation](workflows/solver-validation.yml) | Build optimized and strict-debug, run `make check`, restart and source-mass/PCCU checks, plus a mandatory short Etna two-layer end-to-end test. Optional inputs enable extended Etna restart/1/2/4-thread checks and the independent complete-solver OpenMP comparison. | After source changes and before integrating new numerical developments. |
| [Cell physics and numerical regressions](workflows/physics-regressions.yml) | Run current focused `TESTS/` suites in three groups: state, physical sources and spatial discretization. Each suite has its own saved log. | After modifying a physical law, state conversion or spatial operator. |
| [Static source and documentation audit](workflows/fortran-dependency-audit.yml) | Check Fortran module/Makefile dependencies, Doxygen/formatting, source interface contracts and indexing. Does not compile the full solver. | After modifying module dependencies, APIs, source terms or documentation. |
| [OpenMP performance benchmark](workflows/lie-openmp-benchmark.yml) | Timed end-to-end performance measurements and CPU-resource reports for manually selected grid, horizon and repetition counts. | When assessing OpenMP optimization/scaling, not for every change. |
| [Publish documentation](workflows/publish-documentation.yml) | Test example catalogue generation, build MkDocs + Doxygen (including the API main page) and deploy one GitHub Pages artifact. Only the selected `master` revision may publish. | After merging validated documentation changes into `master`. |
| [Publish Docker image to GHCR](workflows/docker_create_publish.yml) | Build the multi-stage runtime image, check shared-library linkage, and publish `latest` plus an immutable `sha-<commit>` tag with provenance attestation. Only `master` may publish. | After validating and merging changes into `master`; verify package visibility. |

The documentation and Docker publishing workflows can also be run as **build-only previews** on a PR branch. Their deployment and registry push steps are restricted to `master`.

All workflows are intentionally independent. Running the build-only
workflow does **not** trigger any regression, benchmark or publication.

## Cleanup decisions

The historical workflow files for the Intel ifx experiment, D19
settling bisect, two CF gprof investigations, the D20-specific
OpenMP validation and the stateNEW-specific validation were removed
from the active workflow catalogue. Their relevant contemporary
coverage is either handled above or intentionally omitted because
it depended on retired implementation paths.

Unreferenced D17 programs, the obsolete state-equivalence executable and
shell drivers requiring the removed IMEX-stage observer have been removed.
Historical N7 contracts and checkers remain accessible at
[the pre-consolidation revision](https://github.com/demichie/IMEX_SfloW2D_dev/tree/6650582fea44ca5bcbeeb73fc5f631ba0e7bd166/TESTS);
completed evidence and dated acceptance reports remain in `docs/acceptance/`.
The retired packages are not executable coverage of the current source tree.
The active acceptance inventory is validated by
`TESTS/ACCEPTANCE/test_plan_inventory.py`.

Workflows have deliberately narrow permissions: publication workflows
receive the write permissions necessary for their destinations;
other jobs only need read access to the repository.
