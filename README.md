# IMEX-SfloW2D

**IMEX-SfloW2D** is an open-source, depth-averaged **finite-volume solver for two-dimensional geophysical flows**, written in Fortran. It combines a hydrostatic, path-conservative central-upwind spatial discretization with implicit-explicit time integration, ordered physical-source treatments, and OpenMP parallelism.

The code addresses gravity-driven, potentially multiphase flows over non-flat terrain, including granular avalanches, pyroclastic-density-current models, lahars and lava flows. Depending on the configuration, it evolves mixture composition, temperature, pore-pressure or stochastic variables, and exchanges material with the substrate and surrounding atmosphere. The integrated solver supports **one layer** (`N_LAYERS=1`) or **two dynamically coupled layers** (`N_LAYERS=2`), selected directly in the input file.

**Project resources:** [Documentation website](https://demichie.github.io/IMEX_SfloW2D_dev/) · [Fortran API reference](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/) · [Examples](EXAMPLES/) · [Source code](src/) · [Citation](CITATION.cff) · [License](COPYING)

> **Scientific scope:** The two-layer formulation is available for research use within documented restrictions. Its gravitational characteristic spectrum is not always real or complete, and the CFL condition does not detect loss of hyperbolicity. A successful simulation or regression test is not proof of physical validity for arbitrary density contrasts, interlayer shear or constitutive parameters. Read [Two-layer hyperbolicity](docs/TWO_LAYER_HYPERBOLICITY.md) before using coupled results quantitatively.

## Contents

- [Physical model](#physical-model)
- [Physical processes and constitutive laws](#physical-processes-and-constitutive-laws)
- [Numerical methods](#numerical-methods)
- [Code organization](#code-organization)
- [Building from source](#building-from-source)
- [Configuring and running a simulation](#configuring-and-running-a-simulation)
- [Examples](#examples)
- [Output, visualization and restart](#output-visualization-and-restart)
- [Docker](#docker)
- [Verification and reproducibility](#verification-and-reproducibility)
- [Documentation, contributions and citation](#documentation-contributions-and-citation)

## Physical model

### Depth-averaged mixture equations

The solver evolves depth-integrated mass, horizontal momentum, thermal content and constituent inventories on a Cartesian computational grid. Transport and hydrostatic/nonconservative pressure terms are combined with physical source terms for gravity, friction, thermal evolution and mass exchange.

For each present layer `l`, the leading conservative variables per unit horizontal area are

$$
\mathbf q_l =
\left(
M_l,\ M_l u_l,\ M_l v_l,\ M_l c_{p,m,l}T_l,\ \ldots
\right)^{\mathsf T},
\qquad M_l=\rho_{m,l}h_l .
$$

Here `h` denotes layer thickness, `rho_m` mixture density, `u` and `v` depth-averaged horizontal velocities, `T` temperature and `c_p,m` mixture specific heat. **The fourth component is sensible thermal content per unit area, not total mechanical energy.** Additional conservative components represent individual solid and additional-gas species masses and configured transported variables.

Primitive variables contain thickness, velocities, temperature and constituent mass fractions. Conversion between conservative and primitive states uses the configured mixture thermodynamics and admissibility constraints. The state implementation handles wet/dry cells explicitly, and physical laws evaluate cell-local states without owning a second simulation-wide state.

- [State representation, thermodynamics and admissibility](src/state/README.md)
- [Physical equation terms](src/physics/README.md)

### One- and two-layer formulations

| Property | Single layer | Coupled two layers |
| --- | --- | --- |
| Selection | `N_LAYERS=1` | `N_LAYERS=2` |
| Evolving flow | One mixture layer | Lower `L1` and upper `L2` with independent states |
| Shared quantities | Grid, terrain and substrate | Grid, terrain, material/species catalogue and substrate |
| Physical models | Configured for the one layer | Configurable per layer, with pair-wide interlayer laws |
| Transport | Single-layer finite-volume operator | Spatial operator includes interlayer pressure coupling |
| Exchange | Flow–substrate and flow–atmosphere processes | Those external processes plus conservative interlayer mass transfer and momentum drag |

For two layers, the reconstructed surfaces are `eta1 = B + h1` and `eta2 = B + h1 + h2`, where `B` is terrain elevation. The two-layer model is **one coupled system**, not two independent single-layer calculations.

The layer contact/exposure logic distinguishes four local states: both layers wet, only L1 wet, only L2 wet, or both dry. This determines contact with the substrate or atmosphere and the applicability of interlayer processes. Interlayer settling/elutriation transfers material between the layers, together with associated momentum and sensible thermal content, using matched conservative contributions. An independently configured interfacial drag law exchanges horizontal momentum.

The integrated input path needs **no special environment activation flag**. Both modes use the same solver infrastructure, with the scientifically important distinction that the two-layer characteristic problem is only conditionally hyperbolic.

See the [two-layer architecture and governing contracts](docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md), [interlayer physical closures](docs/TL13B_C_PHYSICAL_CLOSURES.md) and [characteristic analysis](docs/TWO_LAYER_HYPERBOLICITY.md).

## Physical processes and constitutive laws

The physical model is assembled from configurable cell-local laws. Availability of a particular numerical treatment depends on the selected source; **not every law supports every treatment**.

| Model family | Implemented capabilities |
| --- | --- |
| Rheology and basal resistance | Voellmy; plastic/yield; lava; lahar; empirical; quadratic; Froude-dependent; Lucas; granular `mu(I)`; combined `mu(I)`–Voellmy laws |
| Material composition | Mixture/carrier properties and multiple solid or additional-gas species, with conservative component transport |
| External mass exchange | Substrate erosion, particle settling/deposition, ambient-air entrainment, carrier-liquid loss, vaporization, and lava deposition, when configured |
| Thermal physics | Temperature-dependent constitutive properties and radiative/convective cooling |
| Pore-pressure physics | Configurable Darcy-related pore-pressure evolution |
| Injection and sources | Configured material-injection mechanisms and boundary/lateral source terms |
| Interlayer mass exchange | L2-to-L1 settling and L1-to-L2 elutriation with carrier co-transfer, donor-inventory and packing limits |
| Interlayer momentum exchange | Pair-wide Chézy-type drag with an analytic relative-velocity update |
| Additional state | Optional transported stochastic and pore-pressure components; substrate and deposition inventories |

Mass exchange with the terrain or atmosphere is **external** to the layer pair; mass exchanged between L1 and L2 is **internal** and must balance across the pair. Physical source groups use shared inventory and admissibility constraints where transfers compete for the same material.

The exact equations, method availability and parameter semantics are documented in the source-family guides:

- [Rheological laws](src/physics/rheology_laws/README.md) and [thermal laws](src/physics/thermal_laws/README.md)
- [Substrate/atmosphere mass exchange](src/physics/mass_exchange/README.md) and [pore pressure](src/physics/pore_pressure/README.md)
- [Interlayer mass exchange](src/physics/interlayer_mass_exchange/README.md) and [interlayer momentum drag](src/physics/interlayer_momentum_exchange/README.md)

A configurable constitutive law is not automatically calibrated for a particular geological event. Use case-specific reference data and parameter validation when interpreting predictions.

## Numerical methods

### Spatial discretization

IMEX-SfloW2D uses a finite-volume formulation on a structured Cartesian grid with terrain geometry, reconstructed face states and central-upwind numerical fluxes. The current **KT / hydrostatic-path PCCU** framework incorporates hydrostatic and nonconservative pressure terms, slope/topographic contributions and wet/dry treatments. For two layers, the same spatial operator also assembles the coupled pressure terms.

Reconstruction, composition normalization, component flux limiting and positivity/admissibility checks are treated explicitly in their respective operators. Central-upwind estimates do not require constructing a complete eigenbasis at every face. **That numerical property is not a proof that the two-layer continuous equations remain hyperbolic.**

### Time integration and physical source terms

Spatial transport is handled explicitly, while physical-source contributions are organized into validated, **ordered Lie groups**. Supported source treatments include:

- **Explicit** source evaluation;
- **Semi-implicit** source evaluation where implemented;
- **Newton/implicit** nonlinear solves, with complex-step Jacobian evaluation and LAPACK/BLAS linear algebra;
- **Analytic** substeps for laws admitting a supported closed-form update.

A treatment name is validated against the selected physical law; unimplemented methods are not silently substituted. REAL and COMPLEX evaluations share physical kernels where required for consistent Newton residuals and Jacobians. Only an admissible full-step candidate is published as the new simulation state.

### CFL time step and parallel execution

The `NUMERIC_PARAMETERS` namelist supports two CFL policies:

- `CFL_POLICY='directional_min'` (default): the minimum of the active directional time-step limits and `MAX_DT`;
- `CFL_POLICY='combined'`: a bound using the sum of active directional wave-speed contributions divided by their grid spacings, also capped by `MAX_DT`.

The two-layer total-depth characteristic-speed envelope is justified only under the conditions described in the [hyperbolicity analysis](docs/TWO_LAYER_HYPERBOLICITY.md). Decreasing the timestep does not repair a non-hyperbolic physical state.

OpenMP parallelism is applied to numerical work over active cells and interfaces, with shared active-domain worksets and reusable thread-local physical workspaces. A representative configuration can set `OMP_NUM_THREADS`, for example `OMP_NUM_THREADS=4`, subject to available hardware and the case's own settings.

For implementation details, see [spatial methods](src/spatial/README.md), [numerical integration](src/numerics/README.md) and the [numerical overview](docs/website/model/numerics.md).

## Code organization

| Directory or file | Responsibility |
| --- | --- |
| [`src/IMEX_SfloW2D.f90`](src/IMEX_SfloW2D.f90) | Main executable, run initialization, advancement and output lifecycle |
| [`src/core/`](src/core/README.md) | Solver-wide configuration, simulation context, model layout and diagnostics |
| [`src/state/`](src/state/README.md) | Canonical conservative/primitive states, conversions, mixture properties and dry-state contracts |
| [`src/geometry/`](src/geometry/README.md) | Grid, topography, slopes and active-domain worksets |
| [`src/spatial/`](src/spatial/README.md) | Reconstruction, numerical fluxes, pressure terms and spatial residuals |
| [`src/physics/`](src/physics/README.md) | Cell-local constitutive laws, physical source rates and exchanges |
| [`src/numerics/`](src/numerics/README.md) | Time stepping, Lie source groups, nonlinear/analytic methods and source validation |
| [`src/io/`](src/io/README.md) | Input parsing, output files, NetCDF and binary restart |
| [`src/stochastic/`](src/stochastic/README.md) | Configured stochastic forcing and reproducible random-state handling |
| [`EXAMPLES/`](EXAMPLES/) | Runnable cases, input files, preparation scripts and case-specific instructions |
| [`TESTS/`](TESTS/) | Unit, numerical, integration and acceptance-oriented checks |
| [`docs/`](docs/README.md) | Scientific architecture, formats, validation, developer guides and website sources |
| [`.github/workflows/`](.github/ACTIONS.md) | Manually dispatched GitHub Actions for verification and publication |

The governing architecture keeps **one authoritative simulation-wide state** in `q/qp(component,layer,x,y)`. Per-layer physical configurations have defined owners; cell-local physical states are reusable workspaces rather than independent persistent copies. Physical kernels do not traverse the grid themselves. See the [development roadmap](docs/IMEX_SfloW2D_Development_Roadmap.md) for current integration status and remaining scientific work.

## Building from source

### Dependencies

A Fortran compiler compatible with the linked libraries, Autoconf, Automake, Make, **NetCDF-Fortran**, and **LAPACK/BLAS** are required. The supplied GNU build is configured for gfortran and OpenMP. Python 3 and case-specific Python libraries may also be needed by example-preparation scripts; they are not required for the Fortran solver itself.

On Ubuntu 24.04, install the standard build dependencies:

~~~bash
sudo apt-get update
sudo apt-get install -y \
    gfortran autoconf automake make \
    libnetcdff-dev libnetcdf-dev liblapack-dev libblas-dev
~~~

Clone and build from the repository root:

~~~bash
git clone https://github.com/demichie/IMEX_SfloW2D_dev.git
cd IMEX_SfloW2D_dev
autoreconf -fi
./configure --with-netcdf=/usr FC=gfortran
make -j2
~~~

The executable is **`src/IMEX_SfloW2D`**. The Autotools-generated configuration files are not committed to the repository.

To use stricter runtime checks in a separate clean build, configure with:

~~~bash
./configure --with-netcdf=/usr FC=gfortran --enable-strict-debug
make clean
make -j2
~~~

This profile enables gfortran diagnostics, bounds checks and floating-point traps. On systems with nonstandard library installations, replace `/usr` with a compiler-compatible NetCDF prefix; `NETCDF_LIBS`, `LAPACK_LIBS` and `OPENMP_FCFLAGS` are configurable build variables.

Run the repository's available local tests with `make check` once the executable is built. See [installation details](docs/website/getting-started/installation.md) and [test workflow descriptions](.github/ACTIONS.md).

## Configuring and running a simulation

### Input-file workflow

The solver reads **`IMEX_SfloW2D.inp` from the current working directory**. Each example includes a suitable input file and documents any DEM, initial-state or Python preprocessing required before execution. Work from the example directory so that its relative paths resolve correctly.

Important input groups include:

| Input group | Purpose |
| --- | --- |
| `RUN_PARAMETERS` | Run identifier, time horizon, output intervals, output options and restart/initialization control |
| `NEWRUN_PARAMETERS` | Grid, layer count (`N_LAYERS`), material classes and physical-model switches |
| `NUMERIC_PARAMETERS` | CFL policy, timestep limits and numerical controls |
| `RESTART_PARAMETERS` / `RESTART_L1_PARAMETERS` / `RESTART_L2_PARAMETERS` | Input raster or conservative-state files for initializing material |
| `INITIAL_PARAMETERS` / `INITIAL_L1_PARAMETERS` / `INITIAL_L2_PARAMETERS` | Initial primitive composition, temperature and velocity |
| `RHEOLOGY_PARAMETERS` / `RHEOLOGY_L1_PARAMETERS` / `RHEOLOGY_L2_PARAMETERS` | Constitutive-law selection, coefficients and source treatment |
| `MASS_EXCHANGE_PARAMETERS` / `MASS_EXCHANGE_L1_PARAMETERS` / `MASS_EXCHANGE_L2_PARAMETERS` | External source and exchange laws |
| `INTERLAYER_PARAMETERS` | Pair-wide interlayer mass and momentum exchange for coupled runs |
| Boundary-condition, thermal and substrate groups | Layer boundaries, thermal evolution and shared terrain inventories |

**Single-layer unsuffixed NAMELISTs are intentional supported input names**, not obsolete compatibility aliases. Two-layer cases use separate L1/L2 groups where applicable. When `RESTART=T`, the listed `.asc` rasters or `.q_2d` snapshots may specify **initial conditions**; they are not necessarily a continuation from a previous simulation. Binary `restart.bin` continuation is a distinct, complete-state mechanism. Two-layer primitive initialization does not currently support nonzero, independent layer release times.

For the authoritative file-format and layer rules, see [initialization and restart](docs/RESTART_INPUT.md). For the actual model and source selections, use an existing [example input](EXAMPLES/) rather than constructing an incomplete NAMELIST from scratch.

### A minimal runnable workflow

The two-layer 1D dam-break is a compact starting point. It uses a flat grid with 200 x-cells and one y-cell; its upper layer starts dry, and interlayer transfer is disabled to isolate transport.

~~~bash
cd EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D
python3 prepare.py
../../src/IMEX_SfloW2D
~~~

Run these commands **after building from the repository root**. The preparation script creates the needed input rasters without running the solver. Outputs are written under `OUTPUT/dambreak_two_layer_1d/`. This example verifies the two-layer input and spatial pathway but does **not** test sustained dynamics with two wet interacting layers.

For OpenMP execution, set the thread count before invoking the solver, for example:

~~~bash
OMP_NUM_THREADS=4 ../../src/IMEX_SfloW2D
~~~

The number of threads actually used also depends on the example's `SERIAL_FLAG` and available OpenMP runtime resources.

## Examples

The [`EXAMPLES/` directory](EXAMPLES/) contains numerical benchmarks and geophysical applications. Each runnable case has its own instructions describing physical assumptions, required data, preprocessing, parameters and expected outputs. Some large terrain DEMs are stored once as ZIP files and must be copied/extracted into the target example's `DEM/` directory.

| Example | Description |
| --- | --- |
| [Two-layer 1D dam-break](EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D/README.md) | Flat-terrain wet/dry propagation; lower layer initially wet, upper layer dry |
| [Etna two-layer case](EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/README.md) | Lower granular collapse and upper particle-laden current with elutriation, atmospheric entrainment and interfacial drag |
| [Vesuvius lahar](EXAMPLES/EXAMPLE_LAHAR/README.md) | Single-layer liquid–sediment flow over volcanic topography |
| [Etna lava](EXAMPLES/EXAMPLE_ETNA_LAVA/README.md) | Temperature-dependent lava flow and fissure input |
| [Taal](EXAMPLES/EXAMPLE_TAAL/README.md) | Radially supplied dilute gas–particle current |
| [BW example](EXAMPLES/EXAMPLE_BW/README.md) | One-dimensional flow with thermal and physical-source effects |

See also the [generated Examples catalogue](https://demichie.github.io/IMEX_SfloW2D_dev/examples/) and [DEM preparation policy](EXAMPLES/DEM_LAYOUT.md). Example parameters are illustrative unless calibration and comparison with observations are expressly documented.

## Output, visualization and restart

### Output layout

Each run writes into a dedicated directory, normally **`OUTPUT/<RUN_NAME>/`** relative to the case directory. The usual output structure includes the following, according to enabled options:

~~~text
OUTPUT/<RUN_NAME>/
  fields/         NetCDF time-dependent physical fields
  diagnostics/    runout and maximum-thickness diagnostics
  snapshots/      time-indexed field or conservative snapshots
  terrain/        terrain grids for interpretation
  probes/         optional probe time series
  vulnerability/  optional vulnerability products
  restart/        complete binary checkpoint (restart.bin)
  metadata/       input backup and snapshot metadata
~~~

`OUTPUT_BASE_DIRECTORY` changes the parent directory. Existing output directories are protected by default: `OUTPUT_OVERWRITE` must be explicitly enabled to reuse an existing run directory, and reuse does not automatically delete old files. Prefer a new `RUN_NAME` for a clean dataset. See [output paths and safety rules](docs/website/getting-started/output-files.md).

### Physical NetCDF fields

With `OUTPUT_NETCDF_FLAG=T`, the principal visualization dataset is **`OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc`**. It contains scalar physical variables on `(x,y,time)` axes, including thickness, velocity, temperature, mixture density, species fractions and configured diagnostics.

Single-layer names are unsuffixed, for example `h`, `u` and `T`. Two-layer quantities are named with `_l1` and `_l2` suffixes, such as `h_l1` and `h_l2`. Shared outputs include terrain `b`, total depth `h_total` and upper free surface `eta_top`. Solid **volume fractions** and **mass fractions** are distinct fields.

**NetCDF is intended for visualization and analysis, not full-state restart.** Continuation of the complete model state uses the current binary checkpoint format. For details, consult [NetCDF field definitions](docs/NETCDF_OUTPUT.md) and [initialization/restart contracts](docs/RESTART_INPUT.md).

### Visualization

The outputs can be examined in ParaView, and terrain/georeferenced imagery can be prepared with QGIS. The [visualization guide](docs/website/getting-started/visualization.md) explains how to display surface height, velocity, temperature and deposition and how to check coordinate-system consistency. The solver's Cartesian coordinates do not, by themselves, define a geographic CRS.

## Docker

The repository includes a multi-stage [Dockerfile](Dockerfile) that builds the solver and packages the runtime executable and Linux libraries without the compiler toolchain. Image publication to **GitHub Container Registry (GHCR)** is a manually triggered workflow; **check that a successful publication exists and that the package is public before attempting to pull it**.

After publication:

~~~bash
docker pull ghcr.io/demichie/imex_sflow2d_dev:latest
cd EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D
python3 prepare.py
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --mount "type=bind,src=$PWD,dst=/work" \
    --workdir /work \
    ghcr.io/demichie/imex_sflow2d_dev:latest
~~~

The image does **not** bundle DEMs or example input files; the mounted case directory supplies them and receives the outputs. For reproducible work, prefer the image's immutable `sha-<40-character-commit>` tag to `latest`.

Alternatively, build from this checkout:

~~~bash
docker build -t imex-sflow2d:local .
~~~

See the [Docker execution guide](docs/website/getting-started/docker.md) and [publication procedure](docs/website/development/publishing.md).

## Verification and reproducibility

The repository includes unit and regression coverage for state conversion, rheology and source laws, REAL/COMPLEX consistency, conservative interlayer exchange, spatial operators, wet/dry behavior, CFL policies, input/output, restart and representative solver cases. Relevant optimized and strict/debug builds and OpenMP comparisons complement these checks.

Run local tests after building:

~~~bash
make check
~~~

**GitHub Actions are intentionally manual**: they do not run automatically on every push. The core validation workflows are:

1. **Static source and documentation audit** — source formatting, Doxygen and dependency/interface checks.
2. **Cell physics and numerical regressions** — state, physical-source and spatial tests.
3. **Solver integration and OpenMP validation** — optimized and strict/debug builds, representative end-to-end simulations and OpenMP/restart checks.

Separate manual workflows publish the combined MkDocs/Doxygen website and the GHCR image. A passing result verifies **the exact Git revision tested**; it does not establish calibration or physical validity for every model regime. Consult the [GitHub Actions guide](.github/ACTIONS.md), [verification overview](docs/website/development/validation.md) and [current roadmap](docs/IMEX_SfloW2D_Development_Roadmap.md).

For reproducible scientific results, record the **commit SHA**, compiler and build flags, complete `IMEX_SfloW2D.inp`, input/DEM checksums, relevant physical model coefficients and runtime settings (including OpenMP threads), together with diagnostics and output data.

## Documentation, contributions and citation

- **User documentation:** [website](https://demichie.github.io/IMEX_SfloW2D_dev/), [build/install](docs/website/getting-started/installation.md), [input/restart](docs/RESTART_INPUT.md), [outputs](docs/NETCDF_OUTPUT.md) and [visualization](docs/website/getting-started/visualization.md).
- **Scientific and architecture references:** [model overview](docs/website/model/index.md), [two-layer architecture](docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md), [hyperbolicity limits](docs/TWO_LAYER_HYPERBOLICITY.md), [source closures](docs/TL13B_C_PHYSICAL_CLOSURES.md) and [development roadmap](docs/IMEX_SfloW2D_Development_Roadmap.md).
- **Developers:** [documentation index](docs/README.md), [Fortran style](docs/FORTRAN_STYLE.md), [adding a rheology](docs/ADDING_A_RHEOLOGY_MODEL.md), [contributing guide](CONTRIBUTING.md) and [Doxygen API](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/).

When using IMEX-SfloW2D in a publication, cite the software according to [`CITATION.cff`](CITATION.cff) and also cite the scientific publications relevant to the selected physical and numerical configuration. The software is authored by **Mattia de' Michieli Vitturi**, Istituto Nazionale di Geofisica e Vulcanologia, with repository contributors as recorded by Git.

IMEX-SfloW2D is distributed under the **GNU General Public License, version 3** ([GPL-3.0-only](COPYING)).
