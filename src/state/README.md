# State: canonical representation of the solver state

## Purpose and role in IMEX-SfloW2D

The `state` directory provides the **canonical representation of the evolving
flow state** and the operations needed to convert, reconstruct, and inspect it.
It is the boundary between the numerical solver, which evolves cell values,
and the physical models, which need temperatures, densities, concentrations,
velocities, and other derived properties.

The directory serves two distinct purposes:

1. **Simulation-wide state ownership.** It owns the conservative and primitive
   arrays on the computational grid, substrate inventories, diagnostic maps,
   and the substrate snapshots used for transactional rollback.
2. **Cell-local state evaluation.** It provides reusable REAL and COMPLEX
   physical-state workspaces and routines that calculate the properties of
   **one cell** from its conservative or primitive variables. These
   evaluations do not require grid indices and are suitable for use inside
   externally parallelized loops.

The numerical algorithms, spatial operators, source-term scheduling, Lie
splitting, and constitutive laws remain outside this directory. In particular,
`state` does **not** own the time-stepping algorithm, determine source rates,
or decide whether to commit or roll back a complete time step. It supplies
the data and conversion/evaluation contracts consumed by those components.

## State representations

The same cell can be described at three complementary levels.

### Conservative variables (`q`)

The time integrator evolves a vector of conserved quantities per unit
horizontal area. Its first four components are

\[
q_1=M=\rho_m h,\qquad
q_2=Mu,\qquad
q_3=Mv,\qquad
q_4=M c_{p,m}T.
\]

Here `h` is flow thickness, `rho_m` is mixture density, `u` and `v`
are depth-averaged horizontal velocities, `T` is temperature, and `c_{p,m}`
is mass-specific mixture heat capacity. The fourth component is **thermal
content per unit area**, not the total kinetic-plus-potential energy.
Additional components carry the masses of solid and additional-gas species,
plus optional transported scalars (such as stochastic and pore-pressure
variables).

### Primitive variables (`qp`)

The primitive representation contains thickness `h`, velocity components
`u` and `v`, temperature `T`, solid and additional-gas **mass fractions**,
and the corresponding optional scalar variables. It is the representation
used for many reconstruction and flux calculations.

For a non-dry cell the conversion uses the mixture composition and
thermodynamic relationships, including

\[
x_{s,i}=q_{s,i}/M,\qquad x_{g,i}=q_{g,i}/M,
\qquad T=q_4/(M c_{p,m}).
\]

The strict conversions and the separate normalization of reconstructed
face compositions are deliberately different contracts.

### Cell-local physical state

`physical_state_real` and `physical_state_complex` gather the primitive
variables and derived properties required by the physics: mixture and carrier
density, reciprocal densities, heat capacities, mass/volume fractions,
reduced gravity, Richardson number, horizontal/total speed, vertical velocity
induced by local bed slope, and curvature acceleration.

These structures are **reusable workspaces**, not additional authoritative
copies of the grid state. Arrays for constituent properties are allocated
and validated during setup; the cell evaluation routines reuse the storage.
The COMPLEX version retains perturbations required by complex-step Jacobians.
No cell index `j` or `k` is passed into the physical-state evaluators.

## State contracts and dry cells

The conversion and physical-state APIs distinguish invalid inputs from
successfully evaluated states through an optional `accepted` output where
provided. For non-dry states they enforce the current contracts on
nonnegative mass/thickness, admissible composition, positive heat capacity
and strictly positive thermal content or temperature. Callers must respect
a failed evaluation rather than interpreting an initialized zero result
as a physical rate or state.

The canonical dry classification uses machine epsilon after rejecting
negative mass or thickness:

- Conservative inputs: `0 <= M <= EPSILON(1.0_wp)`.
- Primitive inputs: `0 <= h <= EPSILON(1.0_wp)`.

A dry physical cell has zero thickness, velocity, transported fractions,
and optional scalars, with ambient thermodynamic properties and
`T=T_ambient`. A dry conservative vector can have `q4=0`; strict
positivity applies to **non-dry** states. Both REAL and COMPLEX constructors
and refresh operations follow this convention.

State-evaluation status and the **global end-of-step admissibility/rollback
decision are distinct**. Source-group Newton evaluation can consume
`accepted` to reject a failed trial without a fatal stop; the time integrator
remains responsible for accepting a complete updated solution.

## Directory contents

| File | Responsibility |
| --- | --- |
| `state.f90` | Public facade that re-exports the canonical state modules. |
| `state_fields.f90` | Owns grid-wide `q`/`qp`, substrate inventories, accumulated diagnostics, allocation/finalization, and snapshots for substrate rollback. |
| `state_types.f90` | Defines reusable REAL/COMPLEX cell-local physical-state types and their initialization, reset, and REAL-to-COMPLEX promotion routines. |
| `state_composition.f90` | Implements strict composition closure, separate face-specific normalization, and conservative mass-budget admissibility checks. |
| `state_mixture.f90` | Evaluates mixture/carrier heat capacities, density, volume fractions, reduced gravity, ideal-gas specific volume, and particle descriptors. |
| `state_conversions.f90` | Converts between `q` and `qp` in REAL/COMPLEX arithmetic; also evaluates dynamic pressure from an already synchronized pair. |
| `state_physical.f90` | Builds complete physical states from `q` or `qp`, refreshes thermodynamics/kinematics, and provides specialized density and temperature evaluations. |
| `inc/` | Shared mathematical kernels included where identical formulas must operate on REAL and COMPLEX quantities. |

### Mathematical kernels in `inc/`

| File | Formula or calculation |
| --- | --- |
| `conservative_from_primitive.inc` | Conservative mass, momentum, thermal content, species and transported scalars from primitive inputs. |
| `temperature_from_energy.inc` | Temperature recovered from conservative thermal content. |
| `velocity_from_momentum.inc` | Horizontal velocities with the established near-dry regularization. |
| `mixture_heat_capacity.inc` | Mass-weighted mixture and carrier heat capacities. |
| `density_from_mass.inc` | Mixture/carrier specific volumes and densities from composition. |
| `mixture_from_mass.inc` | Additional phase properties and volume fractions. |
| `gas_specific_volume.inc` | Ideal-gas specific volume. |
| `dry_thermodynamics.inc` | Canonical ambient thermodynamics for a dry physical cell. |
| `kinematics.inc` | Velocity magnitudes, slope-related vertical velocity and curvature acceleration. |
| `reduced_gravity.inc` | Reduced gravitational acceleration. |
| `richardson.inc` | Depth-averaged Richardson number and its limiting treatment. |

The includes contain **substantive mathematical or physical formulas**,
not wrappers whose only purpose is to call other procedures. Simple local
expressions (for example `xs=qs/M`) remain where they are evaluated;
there is no requirement to extract every repeated assignment into an
include.

## Architectural boundaries and usage guidance

- Grid-wide field ownership is confined to `state_fields`; cell-local
  physical-state calculations operate on values supplied by the caller.
- For routine cell evaluation, start from mass fractions in `qp` or from
  the corresponding conservative cell. `evaluate_mixture_from_volume` is
  reserved for initialization/source input specified in **volume fractions**.
- Use the strict conversion functions for ordinary state operations.
  `normalize_face_composition` belongs to the spatial face-reconstruction
  path and has different clipping/normalization rules.
- Initialize reusable REAL/COMPLEX workspaces before entering cell loops.
  Avoid new allocation in hot per-cell calculations.
- Do not introduce back-references from the state modules to global grid
  indices, time-integrator implementation, or physics source scheduling.
- The legacy `state_2d` and `state_conversion_2d` implementations are
  retired; production code should use the canonical modules in this
  directory.

## Coding conventions and verification

The mandatory [Fortran style policy](../../docs/FORTRAN_STYLE.md) governs
module/procedure organization, four-space lexical indentation, control-block
spacing, 120-column code width and continuation alignment, and Doxygen
documentation (100-column prose). Main operations precede local helpers,
REAL/COMPLEX counterparts remain adjacent, and every procedure has explicit
argument/return documentation.

Formatting and documentation passes must not change mathematical kernels,
the order of floating-point operations within a procedure, source APIs,
OpenMP directives, or simulation behavior. The directory-specific
architecture and work history are recorded in
[STATE_ARCHITECTURE_REDESIGN_ROADMAP.md](../../docs/STATE_ARCHITECTURE_REDESIGN_ROADMAP.md).

Useful repository checks:

```bash
python3 tools/docs/check_fortran_style.py src/
python3 tools/docs/check_fortran_style.py --structural-report src/
python3 tools/docs/check_fortran_style.py --line-width-report src/
python3 TESTS/TEST_SOURCE_DOCUMENTATION/check_documentation.py
python3 .github/scripts/check_fortran_dependencies.py
```

Follow static checks with compilation, the state/equivalence regressions
and representative full simulation examples. Where procedures have been
reordered within a module, compare executable statements **per procedure**:
a whole-file line-order comparison would report the intended reordering
as a difference.
