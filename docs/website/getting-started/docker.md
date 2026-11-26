# Running with Docker

A published Docker image lets you run IMEX-SfloW2D without installing
a Fortran compiler or numerical libraries on the host.

## Get the image

```bash
docker pull ghcr.io/demichie/imex_sflow2d_dev:latest
```

Use a published release-version tag (for example `:v1.2.3`)
or a full-commit `:sha-...` tag when you need reproducible
software. The moving `:latest` tag may change after publication.
Check the [GitHub package page](https://github.com/demichie/IMEX_SfloW2D_dev/pkgs/container/imex_sflow2d_dev)
if the image cannot be pulled.

The image contains the compiled solver and its runtime libraries.
It does **not** contain example directories, DEM rasters or input files.

## Prepare and run a simulation

Prepare the chosen case first, following the
[example instructions](../examples/index.md). From its directory,
which contains `IMEX_SfloW2D.inp` and all required terrain inputs,
run:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/work" \
  --workdir /work \
  ghcr.io/demichie/imex_sflow2d_dev:latest
```

The solver reads the case input from `/work` and writes the
usual `OUTPUT/<RUN_NAME>/` directory into the mounted host
folder. The `--user` option avoids root-owned output on
Linux systems that support `id`.

To use a fixed version, replace `:latest` in the command
with your chosen release or commit tag.

## Build an image locally

From the repository root:

```bash
docker build --tag imex-sflow2d:local .
```

Then use `imex-sflow2d:local` instead of the published GHCR image
in the run command.

After a run, see [Output files](output-files.md) and
[Visualizing results](visualization.md).
