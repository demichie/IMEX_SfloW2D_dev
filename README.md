# IMEX-SfloW2D

**A depth-averaged, finite-volume Fortran solver for volcanic and other
geophysical flows.** The model combines high-resolution spatial
discretization with IMEX time integration, analytic/implicit physical
source treatments and OpenMP parallelism.

The integrated solver supports both **single-layer** (`N_LAYERS=1`) and
**coupled two-layer** (`N_LAYERS=2`) configurations through
`IMEX_SfloW2D.inp`. No special two-layer activation flag is required.
Single-layer carrier configurations include gas, liquid and lava, with
appropriate rheological and mass-exchange laws.

> **Scientific limitation:** the two-layer gravitational characteristic
> spectrum is not unconditionally real. Inverted density ordering and
> some interlayer shear states can produce complex eigenvalues, and
> repeated real eigenvalues need not yield a complete eigenbasis.
> The CFL envelope is not a general hyperbolicity detector. See
> [two-layer hyperbolicity](docs/TWO_LAYER_HYPERBOLICITY.md) and
> [the architecture specification](docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).

## Getting started

### Build from source

Requirements: a compatible Fortran compiler, Autoconf, Automake,
NetCDF-Fortran, LAPACK/BLAS and Make. On Ubuntu 24.04:

```bash
sudo apt-get update
sudo apt-get install -y gfortran autoconf automake make \
    libnetcdff-dev libnetcdf-dev liblapack-dev libblas-dev
```

From the repository root:

```bash
autoreconf -fi
./configure --with-netcdf=/usr FC=gfortran
make -j2
make check
```

The executable is `src/IMEX_SfloW2D`. To build with runtime bounds
checks and floating-point traps, use `./configure --with-netcdf=/usr
FC=gfortran --enable-strict-debug` before `make`.

For custom NetCDF-Fortran installations, replace `/usr` with the
appropriate prefix. The generated Autotools files are intentionally
not tracked in Git.

### Run a prepared example

The one-dimensional two-layer dam-break example builds its input rasters
without running the Fortran solver:

```bash
cd EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D
python3 prepare.py
../../src/IMEX_SfloW2D
```

The solver reads `IMEX_SfloW2D.inp` from the current directory.
Other examples may need a shared DEM ZIP copied and extracted first;
follow each example's own README.

### Run with Docker

Once the [manual GHCR publishing workflow](.github/workflows/docker_create_publish.yml)
has successfully published the container and its package is accessible:

```bash
docker pull ghcr.io/demichie/imex_sflow2d_dev:latest
cd EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D
python3 prepare.py
docker run --rm --user "$(id -u):$(id -g)" \
    --mount "type=bind,src=$PWD,dst=/work" --workdir /work \
    ghcr.io/demichie/imex_sflow2d_dev:latest
```

The Docker image contains the **compiled solver**, not the examples or
data files, and writes results to the mounted case directory. For
reproducibility, use the published `sha-<40-character-commit>` image tag
rather than `latest`. To build the image locally, run
`docker build -t imex-sflow2d:local .` at the repository root.

See the [Docker guide](docs/website/getting-started/docker.md) for
mounts and container permissions.

## Model configuration

- **One layer:** use `N_LAYERS=1`. Unsuffixed physical NAMELIST groups
  are intentional single-layer input forms, not obsolete aliases.
- **Two layers:** use `N_LAYERS=2` and independent L1/L2 configuration
  groups. Both layers share the mesh, terrain and species catalogue.
  Interlayer exchange and drag have one pair-wide `INTERLAYER_PARAMETERS`
  configuration.
- **State:** the solver owns one canonical conservative/primitive
  representation `q/qp(component,layer,x,y)`. Per-layer constitutive
  laws and the accepted-timestep execution are not duplicated.

For input details, see
[initial conditions and restart](docs/RESTART_INPUT.md),
[model architecture](docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md) and
[the model website](https://demichie.github.io/IMEX_SfloW2D_dev/model/).

## Results and restart

By default, generated files are placed under `OUTPUT/<RUN_NAME>/` in
the case directory (override with `OUTPUT_BASE_DIRECTORY`). Principal
subdirectories include `fields/` (NetCDF), `diagnostics/`,
`snapshots/`, `terrain/`, `metadata/` and `restart/`.

- [NetCDF variables and two-layer field names](docs/NETCDF_OUTPUT.md).
- [Binary restart and initial-state inputs](docs/RESTART_INPUT.md).
- [ParaView and QGIS visualization](docs/website/getting-started/visualization.md).

The binary checkpoint has one current development schema; obsolete
checkpoint formats are not maintained.

## Examples

- [Etna two-layer granular collapse and upper dilute current](EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/README.md).
- [Two-layer one-dimensional dam break](EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D/README.md).
- [Vesuvius lahar](EXAMPLES/EXAMPLE_LAHAR/README.md).
- [Complete example catalogue](EXAMPLES/) and its
  [published searchable index](https://demichie.github.io/IMEX_SfloW2D_dev/examples/).

Example parameters are illustrative unless an example explicitly
documents calibration and comparison with reference data.

## Documentation, validation and publication

- [Documentation website](https://demichie.github.io/IMEX_SfloW2D_dev/).
- [Fortran/Doxygen reference](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/).
- [Developer documentation](docs/README.md) and [Fortran style](docs/FORTRAN_STYLE.md).
- [Current tests and scientific acceptance scope](.github/README.md).
- [Contributing](CONTRIBUTING.md).

GitHub Actions are deliberately **manual**: source/documentation audit,
cell-physics regressions, solver/OpenMP integration, documentation
publication and GHCR image publication all run only when explicitly
dispatched. For a new revision, validate the first three workflows
before updating the public site or container. A green build is not a
general scientific validation for all physical states.

## Citation and license

Please cite the software using [CITATION.cff](CITATION.cff) and relevant
scientific literature for your chosen model configuration.
IMEX-SfloW2D is distributed under [GNU GPL v3](COPYING).
