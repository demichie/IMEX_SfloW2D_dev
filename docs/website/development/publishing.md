# Publish documentation and container images

Publication is **intentional and manual**. The workflows do not run on ordinary commits, tags or PRs.

## Documentation website

For a pre-merge check, run **Publish documentation** manually on the PR branch:
MkDocs, Doxygen and example catalogue tests run, but no Pages deployment occurs.
Only the `master` selection uploads and deploys the public site.

1. Merge the validated documentation/source changes into `master`.
2. From **Actions → Publish documentation (MkDocs + Doxygen) → Run workflow**, select `master`.
3. Confirm the **build** and **deploy** jobs succeed.
4. Open the [MkDocs website](https://demichie.github.io/IMEX_SfloW2D_dev/) and [Doxygen API reference](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/).

The build discovers primary simulation guides in `EXAMPLES/`, generates a catalogue without including helper DEM/template README files, and builds Doxygen against the same source SHA. The workflow uploads generation logs and deploys one combined Pages artifact.

## Docker / GitHub Container Registry

For a pre-merge check, run **Publish Docker image to GHCR** manually on the
PR branch: the workflow builds, loads and checks the runtime image locally,
without authenticating to GHCR or publishing tags. Only `master` pushes and attests.

1. Merge the validated Dockerfile/workflow into `master`.
2. From **Actions → Publish Docker image to GHCR → Run workflow**, select `master`.
3. Confirm the image build, runtime-library check, registry push and provenance-attestation steps succeed.
4. Verify the GitHub Packages page for `ghcr.io/demichie/imex_sflow2d_dev`. The repository owner may need to configure package visibility as **public** before anonymous users can pull it.
5. Pull `:latest` for the current release or use its immutable `sha-<commit>` tag for reproducible simulations.

The workflow authenticates using the scoped Actions `GITHUB_TOKEN` with `packages: write`. No personal access token is required in the repository. The image runs the solver against a mounted case directory; it does not run the regression suite by default.

See [Docker execution](../getting-started/docker.md) for pull, mount and run commands.

## Validation before publication

Run the three active source/physics/integration Actions against the exact revision first. A successful build and push do not establish scientific validity of every physical configuration.
