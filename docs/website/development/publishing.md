# Publishing documentation, containers and releases

The project supports **two activation paths**: manual publication from the
selected `master` revision, and automatic, gated publication when a **stable
GitHub Release is published**. An ordinary push, pull request or tag creation
does not publish artifacts.

## Stable GitHub Release

The [release pipeline](https://github.com/demichie/IMEX_SfloW2D_dev/actions/workflows/release.yml)
runs automatically for the `release: published` event, provided that:

- The release is **not a prerelease**.
- Its tag has a version-like name, for example `v1.2.3`.
- The tag resolves to the **current HEAD of `master`**. Older commits
  cannot automatically replace the public `latest` image or website.

Once the gate succeeds, the release pipeline checks **the tagged revision**:
source/style/Doxygen contracts, cell physics, spatial regressions,
optimized/strict-debug solver integration, complete-solver OpenMP comparisons,
and extended two-layer Etna validation. Only after these succeed does it
generate the versioned source/Examples archives and SHA-256 manifest and
attach them to the GitHub Release.

After archive validation, two jobs publish in parallel:
the combined MkDocs/Doxygen website to GitHub Pages, and the tested Docker
runtime to GHCR. The container is tagged `latest`, `sha-<full-commit>` and
with the release version, and receives a provenance attestation.

**A GitHub Release is already visible while its workflow is running.**
If validation or publication fails, the published release is not
automatically deleted; inspect Actions before advertising it as validated.

## Manual documentation publication

The [documentation workflow](https://github.com/demichie/IMEX_SfloW2D_dev/actions/workflows/publish-documentation.yml)
remains manually dispatchable. On the `master` branch it deploys Pages;
on any other selected branch or tag it builds only a preview.

1. In **Actions → Publish documentation (MkDocs + Doxygen)**, select
   **Run workflow** and `master`.
2. Confirm that the example-catalogue tests and `mkdocs build --strict`
   complete successfully.
3. Confirm that Doxygen's **documentation warnings are fatal**, including
   missing/inconsistent parameter descriptions, then inspect its log.
4. Confirm the Pages deployment and open
   [MkDocs](https://demichie.github.io/IMEX_SfloW2D_dev/) and
   [Doxygen](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/).

The Fortran `@mainpage` in `src/IMEX_SfloW2D.f90` is the single Doxygen
main page; `docs/doxygen/mainpage.md` is a supplementary API guide.
Documented private symbols are extracted. The separate source checker
enforces coverage of modules, procedures, arguments and derived types.

## Manual Docker publication

The [GHCR workflow](https://github.com/demichie/IMEX_SfloW2D_dev/actions/workflows/docker_create_publish.yml)
remains manually dispatchable. On non-`master` selections it builds and
**executes a solver fixture inside the local image** without pushing.
On `master` it first smoke-tests the local image, then pushes a matching
build and confirms its registry digest and library linkage.

1. Select **Actions → Publish Docker image to GHCR → Run workflow → master**.
2. Check the TEST_2D runtime smoke test, GHCR digest verification
   and provenance attestation.
3. Confirm package visibility for `ghcr.io/demichie/imex_sflow2d_dev`.
   Public visibility may need an explicit package setting.
4. Pull a revision-specific `sha-...` tag for reproducible results;
   `latest` is a moving tag.

The registry workflow uses the scoped `GITHUB_TOKEN`, not a personal
access token. The image contains the solver and runtime libraries, but no
simulation examples or DEM data; mount the case directory when running.

## Validation and provenance

Manual workflows do **not** implicitly run the release gate. Run the
three independent source, physics and solver Actions on the exact SHA
before a manual scientific publication. A successful CI run validates only
its recorded source revision and test scope, not all two-layer physical
regimes. Source and Examples archives are produced from the **same tag**,
with an inventory and checksums for reproducibility.

See [Docker execution](../getting-started/docker.md) and the
[full workflow catalogue](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/.github/ACTIONS.md).
