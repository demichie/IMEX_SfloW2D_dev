# Running IMEX-SfloW2D with Docker

The GHCR package can be published manually from `master`, or automatically after the validation pipeline for a **published stable GitHub Release**. Both paths build and smoke-test the solver image; check that publication succeeded and that the package is public before pulling.

## Pull the published image

```bash
docker pull ghcr.io/demichie/imex_sflow2d_dev:latest
```

The image contains the compiled solver and Linux runtime libraries. It does **not** contain simulation examples, DEM rasters or input files. Prepare an example on the host first.

## Run a simulation

From an example directory containing `IMEX_SfloW2D.inp` and its required DEM files:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/work" \
  --workdir /work \
  ghcr.io/demichie/imex_sflow2d_dev:latest
```

The executable reads the input from `/work` and writes its normal `OUTPUT/<RUN_NAME>/` tree to the mounted case directory. Supplying the host user/group avoids root-owned output files on Linux and macOS. Adjust the command for platforms that do not provide `id` or a POSIX shell.

The image also receives an immutable revision tag of the form `sha-<40-character Git commit>`. A stable release additionally publishes its version tag (for example `v1.2.3`). Use a revision or version tag instead of `latest` for reproducible results.

## Build locally

From the repository root:

```bash
docker build --tag imex-sflow2d:local .
```

Then replace the GHCR image name in the run command with `imex-sflow2d:local`.

Image build and publication details, permissions and release checks are described under [Publishing](../development/publishing.md). Source-code installation remains available in [Build from source](installation.md).
