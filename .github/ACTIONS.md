# GitHub Actions: manual checks and release automation

Ordinary pushes, pull requests and tag creation **do not** trigger Actions.
Existing workflows remain independently available through **Actions → Run
workflow** (`workflow_dispatch`), including build-only and performance runs.

Publishing a **stable GitHub Release** starts
[Validate and publish a release](workflows/release.yml) automatically
(`release: published`). This workflow uses the existing validation and
publication workflows through `workflow_call`, rather than maintaining
duplicate test implementations.

## Active workflows

| Workflow | Purpose | Activation |
| --- | --- | --- |
| [Build (optimized and debug)](workflows/build-only.yml) | Compile selected optimized or strict-debug profile; no tests | Manual |
| [Static source and documentation audit](workflows/fortran-dependency-audit.yml) | Fortran dependencies, source Doxygen/style checks, acceptance contracts and release-archive inventory | Manual; reusable on release |
| [Cell physics and numerical regressions](workflows/physics-regressions.yml) | State, source and spatial regressions with separate logs | Manual; reusable on release |
| [Solver integration and OpenMP validation](workflows/solver-validation.yml) | Optimized/strict-debug builds, `make check`, restarts, source budgets and mandatory two-layer Etna; optional 1/2/4-thread and extended two-layer tests | Manual; extended checks enabled on release |
| [OpenMP performance benchmark](workflows/lie-openmp-benchmark.yml) | Dedicated end-to-end timing and CPU resource measurements | Manual only |
| [Publish documentation (MkDocs + Doxygen)](workflows/publish-documentation.yml) | Strict MkDocs and Doxygen build, Pages deployment and diagnostics | Manual, or after release validations |
| [Publish Docker image to GHCR](workflows/docker_create_publish.yml) | Local container solver smoke test, GHCR publication, digest verification and provenance attestation | Manual, or after release validations |
| [Validate and publish a release](workflows/release.yml) | Verify tag, run all three suites, build and upload source/example archives, then publish Pages and GHCR | Automatically on a published stable release |

## Stable-release gate

A release is eligible for automatic publication only if it is **not a
prerelease** and its version-like tag resolves to **the current `master`
HEAD**. This prevents the `latest` container tag and public documentation
from silently reverting to a previous commit.

1. Publish the GitHub Release for a tag at the currently selected `master`
   commit (for example `v1.2.3`). Creating a tag without publishing a
   Release does not start the pipeline.
2. Audit, physics regressions and optimized/strict-debug solver integration
   run on the **release tag SHA**. Extended two-layer and end-to-end OpenMP
   checks are enabled for stable releases.
3. Only after all validation succeeds, create the versioned source and
   Examples archives with a SHA-256 manifest, and attach them as release
   assets.
4. Only after archive validation succeeds, publish MkDocs/Doxygen and the
   tested container with `latest`, `sha-<full-commit>` and a version tag.
   Docker provenance attestation uses the pushed digest.

A failed gate or regression blocks both publication jobs; an already
published GitHub Release is **not** automatically withdrawn. Review its
status before promoting it as a validated scientific version.

## Manual publishing and previews

The two publishing workflows still support manual dispatch on `master`.
For a different branch or tag, they perform a build-only preview and do
**not** deploy to Pages or push GHCR. Manual publication is an explicit
maintainer decision; unlike the release pipeline it does not invoke CI
automatically. All outcomes belong to their exact tested SHA.

Doxygen is generated with strict warning failure and private procedures
visible. The separate source documentation checker covers modules,
procedures, argument directions, derived types and components.

## Permissions

All validation workflows use read-only repository permissions. The release
orchestrator grants `contents: write` only to the archive-upload job,
`pages: write` and `id-token: write` only to website publication, and
`packages: write`/`attestations: write`/`id-token: write` only to container
publication. No repository secrets or personal tokens are required.

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
