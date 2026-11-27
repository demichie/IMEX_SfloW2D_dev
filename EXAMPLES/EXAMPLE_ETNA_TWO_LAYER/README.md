# Etna two-layer — granular collapse and upper dilute current

This example represents the collapse and downslope motion of a concentrated
granular deposit (lower layer, L1), which feeds an initially dry upper
particle-laden gas current (L2). The upper current can entrain ambient air;
the lower current can deposit particles onto the substrate. Interlayer
mass transfer and interlayer momentum drag are separate processes.

The source-law coefficients and single-size particle population are
illustrative, not a calibrated reconstruction of an Etna event.
The coupled model's characteristic limits, including complex eigenvalues
for some density and velocity states, are explained in
[two-layer hyperbolicity](../../docs/TWO_LAYER_HYPERBOLICITY.md).
The simulation runs directly with `N_LAYERS=2`.

Run every command from **EXAMPLES/EXAMPLE_ETNA_TWO_LAYER**. Preprocessing
scripts and input files are local to this example. The only file shared
with `EXAMPLE_ETNA` is the original DEM ZIP, which is not tracked twice.

## 1. Copy the shared DEM once

From the repository root:

```bash
cd EXAMPLES/EXAMPLE_ETNA_TWO_LAYER
cp ../EXAMPLE_ETNA/DEM/Etna2014_crop.zip DEM/
unzip DEM/Etna2014_crop.zip -d DEM/
```

The DEM extraction is a manual prerequisite, not part of `prepare.py`.
Check that `DEM/Etna2014_crop.asc` exists before proceeding.

Requirements for preprocessing: Python 3, NumPy and pandas. The compiled
solver is expected at `../../src/IMEX_SfloW2D`.

## 2. Prepare topography and initial thickness (no solver run)

```bash
python3 prepare.py
```

`prepare.py` never starts the Fortran program. It uses the local
`create_input_ellipsoid.py` and `input_ellipsoid.py` and produces:

- `DEM/Etna2014_crop.asc`: manually extracted original DEM (input).
- `DEM/new_dem.asc`: terrain modified by the excavation.
- `DEM/init_volume.asc`: initial L1 thickness raster.
- `DEM/ellipsoid.asc`: excavation geometry for inspection.

The ellipsoid uses the Etna center `(500430, 4177500)` and a target
volume of `100000 m³`. When changing the DEM or excavation parameters,
delete generated rasters before rerunning preprocessing. No other
example's preprocessing scripts are needed.

## 3. Run and interpret the two-layer simulation

The solver reads the *primitive* initial conditions directly from
`IMEX_SfloW2D.inp` and uses its normal mixture and conservative-state
conversions. A preliminary single-layer simulation or external
conservative-state initialization is **not** needed.

```bash
../../src/IMEX_SfloW2D
```

### 3.1 Initial state and numerical setup

| Property | Lower layer (L1) | Upper layer (L2) |
| --- | --- | --- |
| Interpretation | Concentrated granular avalanche | Particle-laden gas current forming during motion |
| Initial thickness | `DEM/init_volume.asc` | `NONE` (initially dry) |
| Initial solid volume fraction | 50% | 0% (dry) |
| Initial temperature | 373 K | 373 K specified, but unused while dry |
| Initial horizontal velocity | 0 m/s | 0 m/s |
| Rheology | Voellmy: `MU=0.3`, `XI=300` | Quadratic resistance: `FRICTION_FACTOR=0.02` |

The run uses `N_LAYERS=2`, `GAS_FLAG=T`, `LIQUID_FLAG=F`,
`N_SOLID=1` and `N_ADD_GAS=0`. The single solid species has
`RHO_S=2000 kg/m³` and `DIAM_S=1.0E-3 m` (1 mm). Thus L2 is
*not* a calibrated fine-ash population; both layers share the same
solid species. The source temperatures and ambient properties include
`T_AMBIENT=300 K` and `PRES=101300 Pa`.

The current input sets `T_END=50 s`, `DT_OUTPUT=2 s`, `DT0=0.01 s`,
`MAX_DT=1 s` and `CFL=0.24`, with a `250 x 150` grid of
10 m cells. `KT_TRACE_RECONSTRUCTION='froude_blend'` is selected.
Results are written to `OUTPUT/exampleEtna_two_layer/`.

The independent `RESTART_L1_PARAMETERS` and `RESTART_L2_PARAMETERS`
groups define the initial thickness files, not a preliminary simulation:
`RESTART_FILES` belongs only to L1 and `RESTART_FILES_L2` only to
L2. `'NONE'` denotes an initially dry layer. The accompanying
`INITIAL_L1_PARAMETERS` and `INITIAL_L2_PARAMETERS` specify the
primitive composition, velocity and temperature. Both layers can
alternatively start wet, but nonzero independent release times are not
supported in the coupled model. See [restart input](../../docs/RESTART_INPUT.md).
Continuation from a saved conservative state remains a separate restart
mechanism.

### 3.2 Active mass pathways

The following table describes the **current input**, not every process
supported by the solver.

| Pathway | Current setting | Physical and numerical effect |
| --- | --- | --- |
| **L1 -> substrate** | **ON**: `MASS_EXCHANGE_L1%DEPOSITION%SETTLING_FLAG=T` | Settling removes mobile L1 solids, adds to the substrate deposition inventory, and can withdraw associated carrier according to the deposition/porosity closure. |
| Substrate -> L1 | OFF: `MASS_EXCHANGE_L1%EROSION%COEFFICIENT=0.0` | No prescribed substrate erosion replenishes the flow. |
| Atmosphere -> L1 | OFF: `MASS_EXCHANGE_L1%ENTRAINMENT%FLAG=F` | The lower layer does not directly entrain ambient air. |
| **L1 -> L2** | **ON**: `INTERLAYER_MASS_EXCHANGE%ELUTRIATION_FLAG=T` | Shear-driven elutriation transfers solid mass *and gas already present in L1* to the upper layer. |
| L2 -> L1 | OFF: `INTERLAYER_MASS_EXCHANGE%SETTLING_FLAG=F` | No particle settling back from L2 to L1 is applied. |
| **Atmosphere -> L2** | **ON**: `MASS_EXCHANGE_L2%ENTRAINMENT%FLAG=T` | Ambient-air incorporation adds gas/mixture mass and sensible thermal energy to a wet L2, changing its dilution and thermodynamic state. |
| L2 -> substrate | OFF: `MASS_EXCHANGE_L2%DEPOSITION%SETTLING_FLAG=F` | No direct basal deposition from the upper layer; interlayer settling, if enabled later, would first transfer solids to L1. |
| Substrate -> L2 | OFF: `MASS_EXCHANGE_L2%EROSION%COEFFICIENT=0.0` | No substrate erosion source for L2. |

The contacts are determined per cell from the current layer thicknesses:
L1 contacts the ground whenever wet; L2 is the atmospheric-facing layer
when wet. When L2 is absent, L1 is exposed to the atmosphere, although
its own entrainment flag is disabled here. The L2 air-entrainment law
also requires positive thickness and positive mixture Richardson number.
Enabling a source flag does not imply a nonzero rate in every cell.

### 3.3 L1 deposition and fixed topography

L1 deposition uses the shared particle-settling kernel and the current
L1 composition, rheology, carrier properties and timestep. The removal
is a source of **negative mobile solid mass**, paired with positive
substrate deposition inventory. The carrier withdrawal, if any, follows
the substrate-porosity closure; it is not the same process as transferring
gas from L1 into L2. The mass-exchange coordinator applies inventory
and packing constraints to the local source transaction.

`SUBSTRATE%TOPO_CHANGE_FLAG=F` leaves the topographic surface used by
the flow solver unchanged during the simulation. Deposition can still
be tracked in the substrate inventories: a recorded deposit does not
automatically imply an evolving terrain elevation in this example.
The example does not activate terrain erosion.

### 3.4 L1 -> L2 elutriation and dynamically associated gas

L2 starts dry and is populated by the `LINEAR_USTAR` interlayer
elutriation law. Its current settings are:

```fortran
INTERLAYER_MASS_EXCHANGE%ELUTRIATION_FLAG=T,
INTERLAYER_MASS_EXCHANGE%ELUTRIATION_SOLID_VOLUME_FRACTION=-1,
INTERLAYER_MASS_EXCHANGE%ELUTRIATION_CARRIER_MODE='LIMIT_SOLIDS',
INTERLAYER_MASS_EXCHANGE%ELUTRIATION_LAW='LINEAR_USTAR',
INTERLAYER_MASS_EXCHANGE%ELUTRIATION_EFFICIENCY=0.01,
INTERLAYER_MASS_EXCHANGE%ROUSE_CRITICAL=100.0,
INTERLAYER_MASS_EXCHANGE%SETTLING_FLAG=F,
INTERLAYER_MASS_EXCHANGE%METHOD='EXPLICIT',
INTERLAYER_MASS_EXCHANGE%GROUP=4,
```

The elutriation law proposes solid-volume fluxes using interfacial
shear, solid availability and the Rouse filter. A critical Rouse
number of `100` is intentionally permissive and requires scientific
calibration; it is not a validated value for Etna.

**`ELUTRIATION_SOLID_VOLUME_FRACTION=-1` is a mode selector, not a
negative concentration.** It derives the transferred parcel's solid
volume fraction from the *current local L1 donor inventories*, after
spatial transport and any preceding Lie source updates:

$$
\alpha_{\mathrm{lift}}
  =\frac{V_{s,1}}{V_{s,1}+V_{c,1}},\qquad
E_c=\frac{1-\alpha_{\mathrm{lift}}}{\alpha_{\mathrm{lift}}}E_s.
$$

Here `V_{s,1}` and `V_{c,1}` are local L1 solid and carrier volumes
per unit area; `E_s` is the accepted upward solid-volume rate and
`E_c` is the associated upward carrier-volume rate. With
`LIMIT_SOLIDS`, the requested parcel composition is preserved by
limiting solid transfer when donor gas is insufficient. The paired
limiter also enforces donor inventories and final solid packing.
If L1 has no solids available, there is no upward co-lift.

Specifying a positive fraction such as `0.15` would instead prescribe
the lifted parcel composition at **15% solids and 85% gas by volume**,
regardless of the mean local L1 composition. Omitting the field
entirely is equivalent to `-1`.

Solid and carrier transfers subtract mass from L1 and add the same
mass to L2. They also carry donor-state horizontal momentum and
sensible thermal energy. This **internal** transfer conserves the
L1+L2 sums of solid species, carrier mass, momentum and thermal energy
for the paired source operation; it does not create atmospheric gas.
The resulting L2 concentration is *not* constrained to equal
`alpha_lift`, because L2 also undergoes transport and air entrainment.

### 3.5 Ambient-air entrainment into L2

`MASS_EXCHANGE_L2%ENTRAINMENT%FLAG=T` activates an independent
**external** source, with `METHOD='EXPLICIT'` and `GROUP=4`.
For a wet, eligible L2 with positive Richardson number, the source
adds ambient-air mass at ambient density and thermal energy using
`T_AMBIENT=300 K`. It does **not** directly add solid mass or
horizontal momentum. Consequently, it can dilute and cool L2 and
alter its velocity through the changed mixture mass. This is distinct
from the gas carried upward out of L1 by elutriation.

`N_ADD_GAS=0` means no additional gas species is separately tracked
in this setup. The gas carrier is represented through the mixture
mass and solid-composition variables.

### 3.6 Interlayer momentum drag and source scheduling

The two layer rheologies are independent: L1 uses Voellmy friction
and L2 uses quadratic resistance (`FRICTION_FACTOR=0.02`), both
with Newton treatment in Lie group 1. Thermal sources use group 2.

The **interlayer momentum** source is separately configured as
`MODEL='CHEZY'`, `CHEZY_DRAG_COEFF=0.03`, `METHOD='ANALYTIC'`
and `GROUP=5`. It acts only where both layers are thicker than
1 mm, transferring equal and opposite momentum impulses according
to their relative velocities. It does not transfer additional mass.
The drag dissipates resolved relative kinetic energy **without adding
that lost mechanical energy to either layer's sensible heat**. This is
an acknowledged, accepted limitation of the present physical model:
the drag-only source conserves pair momentum but does not conserve
combined kinetic-plus-thermal energy. No frictional-heating law is
included or implied.
The Chézy interlayer coefficient must not be confused with L2's
quadratic rheological coefficient; neither is calibrated here.

L1 substrate exchange, L2 atmospheric exchange and paired interlayer
mass exchange all use **explicit Lie group 4**. They are distinct
physical sources evaluated within the configured ordered Lie
sequence; group 5 subsequently applies interlayer drag. Advecting
solids and gas across the interface also transports momentum and
energy, independently of this drag-only impulse.

### 3.7 Mass budgets and interpretation of solver output

For the currently enabled *mass sources*, the paired L1 -> L2
elutriation cancels when summing both layers. Omitting boundary
transport and any other optional mass-removal mechanism, the net
external source for the combined moving flow is:

$$
\dot M_{\mathrm{mobile},1+2}
  = \dot M_{\mathrm{air\ into\ L2}}
    - \dot M_{\mathrm{L1\ deposited}}.
$$

The total mobile **solid** inventory decreases through L1
deposition, while solid elutriation merely redistributes solids
between the layers. Gas entrainment from the atmosphere increases
mobile gas mass. Spatial fluxes through the domain boundaries and
other enabled numerical/physical removal paths must be considered
when auditing a full run.

**Important:** the periodic terminal summary labeled `mass`,
`solid mass`, `volume` and `area` currently reports **L1 only**
(`src/IMEX_SfloW2D.f90`); these are *not* whole-column L1+L2
totals. Use layer-resolved output to examine L2 formation and compute
combined mass balances. In particular, falling L1 mass can reflect
both external deposition and internal transfer to L2.

## 4. Optional checks — commands from this directory

Unit checks (no solver or DEM required):

```bash
python3 -m unittest discover -s . -p 'test_*.py'
```

Short controlled integration smoke:

```bash
python3 run_validation.py \
  --solver ../../src/IMEX_SfloW2D \
  --workdir ./validation-smoke
```

Extended conservation/restart/OpenMP checks:

```bash
python3 extended_validation.py \
  --solver ../../src/IMEX_SfloW2D \
  --workdir ./validation-extended
```

The checkers run from this example and create isolated, new or empty
work directories. They copy the already extracted local DEM for preprocessing but **do
not** invoke a preliminary single-layer solver to construct initial
snapshots. The two-layer executable alone constructs its conservative
initial state when each test begins.

The short test checks snapshots, ESRI/NetCDF output, budgets and binary
restart presence. The extended test compares continuous vs restarted
trajectories, 1/2/4 OpenMP threads and L2 formation over a 2 s
trajectory. These checks cover outputs, budgets and reproducibility;
they do not establish physical calibration or guarantee that all
states have real characteristic eigenvalues.

## Source files

`prepare.py` generates only DEM rasters. `IMEX_SfloW2D.inp` defines the
two-layer run, including its primitive initial conditions. The local
`input_ellipsoid.py` and `create_input_ellipsoid.py` generate the
original Etna excavation without using another example's scripts.
`run_validation.py`, `extended_validation.py`, and `test_*.py`
provide optional regression evidence.

The `.asc` rasters and `OUTPUT/` are generated and untracked. The
only copied source data is `DEM/Etna2014_crop.zip`.

## Positivity and Lie-source diagnostics

The timestep preparation computes a spatial mixture-mass positivity
bound `dt_pos` for both layers and combines it with the CFL bound.
It does not independently guarantee positivity of each solid/gas
component or admissibility of explicit Lie mass-exchange sources.
To investigate a rejected source update, enable
`IMEX_LIE_DIAGNOSTICS=1`; timestep selection can be inspected with
`IMEX_TIMESTEP_DIAGNOSTICS=1`. In particular,
`TL07: paired interlayer physical transfer rejected` identifies
a failure in the interlayer mass-source construction, not a direct
failure of the spatial `dt_pos` calculation.

On a rejected spatial update, the solver prints `KP_POSITIVITY` records
after `LIE_SPATIAL_FAILURE` when `N_LAYERS=2`. They report the
cell-owned W/E/S/N transport depths, their four-face mean (equation 3.31
of Kurganov–Petrova), conservative face masses obtained from the same
selected traces as the PCCU mass flux, and the characteristic speed
bounds and CFL values. The diagnostic also reconstructs the eight mass
terms analogous to equation 3.32, distinguishing the four cell-owned
terms from the four incoming neighbor terms. `conservation_offset`
reports the difference between the published cell mass and the mean
of the reconstructed face masses; a nonzero offset means the eight-term
positivity argument does not apply without that extra contribution.

The eight-term decomposition assumes the unmodified central-upwind
mass flux. Dry-donor suppression and a degenerate-speed fallback can
alter it; compare the reported decomposition with the actual
`LIE_SPATIAL_FAILURE` mass update and `LIE_FACE_FAILURE` fluxes.
This is read-only failure evidence, not a positivity fix.

## L2 face-base and reconstruction audit

The failure-only `KP_BASE` diagnostic uses the completed L1 reconstruction as
the effective L2 base. It prints, in W/E/S/N order, each cell-owned base,
the opposite reconstructed base at the same interface, their jump
(`opposite - own`), and the difference between the cell-owned base
average and `B_cent + h1_cell`. These are distinct questions:
a jump between adjacent cells is allowed by the current reconstruction
contract, while the cell-owned four-face average must still reproduce
the cell-center effective base.

The `KP_TRACE` records separate first-pass primitive reconstructions
(`raw_h`, `raw_mass`, `raw_density`) from the selected PCCU transport
traces (`KP_POSITIVITY h_W,E,S,N`, `mass_W,E,S,N`, and
`selected_density`). Both mass sets use the canonical
`primitive_to_conservative` conversion. Compare their differences
from the published conservative cell mass to determine whether a
large reconstruction offset is already present before the L2
hydrostatic/transport selection. The test
`TESTS/TEST_SPATIAL_OPERATOR` also checks the effective base's
cell-owned average identity and independently retained face jumps.

These diagnostics do not prove that the two-layer pressure system
satisfies the Kurganov-Petrova theorem; they test the necessary mass
reconstruction and numerical speed hypotheses without modifying the
numerical method or mass-exchange sources.

## Source integration

The table follows ascending **Group** numbers (the order of Lie source updates after spatial transport). **Scope** identifies the layer or layer pair affected; **Method** is the configured numerical treatment; **Status** distinguishes enabled sources from settings that do not activate a process; **Role** explains the physical contribution. An entry marked *Configured only* is not an executed source group. Group 0 denotes the initial explicit source stage, not a time-integration stage count.

| Group | Scope | Method | Status | Role |
| ---: | --- | --- | --- | --- |
| 0 | L1 + L2 | `EXPLICIT` | Initial stage | No optional curvature or injection source is enabled. |
| 1 | L1 | `NEWTON` | Active | Voellmy basal friction modifies L1 momentum. |
| 1 | L2 | `NEWTON` | Active | Quadratic resistance modifies L2 momentum. |
| 2 | L1 + L2 | `NEWTON` | Configured only | Thermal cooling is not enabled. |
| 4 | L1 | `EXPLICIT` | Active | Solid deposition removes material from L1 into the substrate. |
| 4 | L2 | `EXPLICIT` | Active | Ambient-air entrainment adds material to L2. |
| 4 | Interlayer | `EXPLICIT` | Active | Elutriation transfers solids and associated carrier between layers. |
| 5 | Interlayer | `ANALYTIC` | Active | Chézy interlayer drag exchanges momentum between the layers. |
