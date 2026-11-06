# Physics: cell-local physical models and source terms

## Purpose and role in IMEX-SfloW2D

The `physics` directory implements the **physical laws and physical
contributions to the depth-averaged multiphase equations**. It contains
constitutive rheologies, thermal and pore-pressure models, mass exchange
with the substrate or carrier, mass-injection mechanisms, and the
equation-level flux and pressure-coefficient evaluations used by the
spatial discretization.

Its central responsibility is to **evaluate physics, not to advance the
simulation**. The physical kernels return rates, fluxes, or local
analytical results from the state and geometric information supplied
by their callers. The numerical algorithms decide when and how these
results enter a time step.

In the current architecture:

- [`state`](../state/README.md) owns the canonical conservative and
  primitive representations, reusable REAL/COMPLEX cell-local physical
  states, and the mixture-property calculations needed by physics.
- **`physics`** implements the laws and assembles their cell-local
  contributions, with dispatchers for physical models and supported
  numerical treatments.
- **`numerics`** owns the source catalog, validated Lie groups, Newton
  solver, group execution, and time integration.
- **`spatial`** owns finite-volume reconstruction and spatial
  operators. It calls physical equation-level evaluations where
  needed; the physical modules do not own the computational mesh
  traversal.

This separation keeps the numerical methods independent of the details
of any single constitutive law and supports wider OpenMP loops, future
MPI decomposition, GPU-oriented execution, and alternative physical
evaluators without embedding grid indices in the physical kernels.

## Physical contributions to the equations

The model evolves conservative cell quantities such as mixture mass,
horizontal momenta, thermal content, and transported constituent
masses. For a wet cell, the leading variables are

\[
\mathbf{q} =
\left(M,\; Mu,\; Mv,\; M c_{p,m}T,\;\ldots\right)^\mathsf{T},
\qquad M=\rho_m h.
\]

The fourth component represents **thermal content per unit area**;
it is not kinetic or gravitational energy.

The physical modules cover two related categories:

1. **Equation-level transport physics:** local inertial fluxes,
   hydrostatic coefficients, characteristic-speed information, and
   limits on species transport. These are evaluated in
   `equation_terms_2d.f90` for use by the spatial methods.
2. **Non-transport physical contributions:** rheological momentum
   exchange, heating/cooling, pore-pressure evolution, material
   injection, entrainment, deposition, erosion, vaporization, carrier
   loss and lava deposition. These are evaluated in the corresponding
   physical model modules and coordinated by `source_terms_2d.f90`.

The modules are not a second copy of the state representation. They
consume conservative/primitive vectors or the reusable
`physical_state_real` and `physical_state_complex` objects provided
by the state layer.

## Source-model interfaces

Source models are organized around **four physical evaluation
interfaces**. Each model family retains the relevant entry points,
including explicit rejection where a treatment is not implemented.

| Interface | Purpose |
| --- | --- |
| **Explicit** | Evaluate a REAL-valued source rate from the current local state. |
| **Semi-implicit** | Evaluate the REAL-valued contribution treated semi-implicitly by the time integrator. |
| **Implicit** | Evaluate a source for Newton residuals and the COMPLEX perturbations used in complex-step Jacobians. REAL implicit evaluators are provided where needed by the residual path. |
| **Analytic** | Integrate a supported physical contribution over a prescribed interval, avoiding a Newton solve where an appropriate analytic update exists. |

These are **physical interfaces**, not instructions to execute all four
treatments for every source. Availability and scheduling are validated
by the immutable catalog in
[`numerics/source_catalog_2d.f90`](../numerics/source_catalog_2d.f90)
and the group registry in
[`numerics/source_registry_2d.f90`](../numerics/source_registry_2d.f90).
Unsupported treatments must not silently return a zero contribution
as though the law had been evaluated.

### REAL and COMPLEX evaluation

REAL implementations are used for explicit and semi-implicit
contributions and for the REAL residual of a Newton group. COMPLEX
implementations propagate the imaginary perturbation needed for a
complex-step Jacobian. Where a mathematical contribution is genuinely
identical between the two precisions, a small `.inc` kernel carries
the shared formula; it is not a forwarding-only wrapper.

This distinction is especially important in rheology: explicit rates
may be assembled from established implicit and semi-implicit
contributions where that avoids duplication, but **semi-implicit rates
remain REAL**. The source interfaces and their implementations must
preserve the selected numerical treatment rather than merely matching
a common name.

### Source assembly and numerical execution

`source_terms_2d.f90` is the physical coordinator for source
evaluation. It dispatches the selected laws, assembles conservative
rates, handles coupled source limits, and exposes separate REAL and
COMPLEX implicit paths. Its explicit transaction assembly returns
candidate contributions and status; it does not directly commit the
simulation-wide state.

The broader source workflow is:

1. **Configuration:** `source_configuration_2d.f90` translates
   physical settings into names and validated source-group metadata.
2. **Evaluation:** `source_terms_2d.f90` and the specialized
   submodules evaluate one cell's physical rates or analytical result.
3. **Numerical solution:** the source executor and Newton/analytical
   solvers in `numerics` apply the configured treatments and control
   the group progression, convergence, and rollback.

The global end-of-step admissibility decision belongs to the
time-integration machinery. Local state reconstruction may report
`accepted=.FALSE.`, including during Newton trials; this status must
be handled by the caller rather than turning an invalid evaluation
into a valid zero source. The canonical dry-state convention itself
is defined by the `state` layer.

## Directory contents

### Top-level physical modules

| File | Responsibility |
| --- | --- |
| `equation_terms_2d.f90` | Local inertial transport flux, hydrostatic pressure coefficient, Cartesian characteristic-speed estimates, and component-mass-flux limiting for spatial operators. |
| `source_configuration_2d.f90` | Convert physical configuration flags and model selections into the validated source metadata used by the numerical source registry; no per-cell rate evaluation. |
| `source_terms_2d.f90` | Assemble and dispatch physical source contributions, coupled explicit limitations, analytical source updates, REAL Newton residual rates and COMPLEX Jacobian rates. |

### `rheology_laws/` — momentum sources

Rheological laws calculate momentum exchange and friction-related
updates. The dispatcher `rheology_2d.f90` selects a law by its
canonical identity and exposes the explicit, semi-implicit, implicit
and analytic paths supported by that law. It also provides the
friction-only momentum limiting used to avoid unphysical reversal or
rotation of momentum by a dissipative substep.

| File | Responsibility |
| --- | --- |
| `rheology_2d.f90` | Model dispatch, common rheology entry points, curvature-related helpers and friction-only momentum limiting. |
| `rheology_analytic_utils_2d.f90` | Mathematical utilities for analytic velocity updates, including scalar root solvers, Lambert W and exact frozen-coefficient flows. |
| `rheology_parameters_2d.f90` | Shared input storage for rheological coefficients used by several models. |
| `rheology_voellmy_2d.f90` | Voellmy-type resistance. |
| `rheology_plastic_2d.f90` | Plastic/yield-related resistance. |
| `rheology_lava_2d.f90` | Temperature-dependent lava rheology. |
| `rheology_lahar_2d.f90` | Lahar rheological resistance. |
| `rheology_empirical_2d.f90` | Empirical resistance law. |
| `rheology_quadratic_2d.f90` | Quadratic drag law. |
| `rheology_froude_2d.f90` | Froude-dependent resistance. |
| `rheology_lucas_2d.f90` | Lucas rheology. |
| `rheology_mui_2d.f90` | \(\mu(I)\)-related rheology. |
| `rheology_mui_voellmy_2d.f90` | Combined \(\mu(I)\) and Voellmy formulation. |
| `rheology_newmodel_2d.f90` | Additional selectable rheology model. |

The directory also contains
`rheology_empirical_drag.inc`,
`rheology_lahar_drag.inc`,
`rheology_lava_viscous.inc`,
`rheology_newmodel_linear.inc`,
`rheology_plastic_yield.inc`,
`rheology_quadratic_drag.inc`, and
`rheology_voellmy_turbulent.inc`. These hold shared mathematical
contributions used in REAL/COMPLEX law evaluations.

### `mass_exchange/` — exchange between mobile flow and surroundings

Mass exchange can modify mixture and species inventories, momentum,
thermal content, and substrate fields. Because several processes
compete for the same available material, their contributions cannot
always be limited independently.

| File | Responsibility |
| --- | --- |
| `mass_exchange_2d.f90` | Coordinate cell-local exchange processes, aggregate conservative source rates, and enforce shared carrier-inventory and packing constraints. |
| `mass_exchange_erosion_2d.f90` | Substrate erosion and entrainment of eroded material into the flow. |
| `mass_exchange_deposition_2d.f90` | Deposition of mobile constituents. |
| `mass_exchange_entrainment_2d.f90` | Entrainment exchanges. |
| `mass_exchange_carrier_liquid_loss_2d.f90` | Loss of mobile liquid carrier. |
| `mass_exchange_vaporization_2d.f90` | Vaporization-related mass and energy exchange. |
| `mass_exchange_lava_deposition_2d.f90` | Lava deposition and associated substrate update contributions. |

`mass_exchange_2d.f90` uses reusable, cell-local
`mass_exchange_scratch_type` storage. Its coordinated limiters protect
shared mass inventories and packing conditions. **Rates and proposed
changes are not themselves a global state commit.** Currently
unsupported semi-implicit, implicit or analytic exchange entry points
remain explicit errors rather than silently approximating an
unimplemented treatment.

### `mass_injection_terms/` — imposed sources and boundary injection

This family provides bottom circular, bottom fissural and lateral
radial material injection, and boundary injection traces. It uses
supplied local geometric quantities and precomputed injection-material
properties.

| File | Responsibility |
| --- | --- |
| `mass_injection_2d.f90` | Dispatch and assemble conservative injection contributions across the supported source geometries. |
| `mass_injection_coefficients_2d.f90` | Evaluate physical-time injection-window coefficients and prepare material properties shared across cells. |
| `mass_bottom_circular_injection_2d.f90` | Bottom circular source. |
| `mass_bottom_fissural_injection_2d.f90` | Bottom fissural source. |
| `mass_radial_injection_2d.f90` | Lateral radial source. |
| `mass_boundary_injection_2d.f90` | Primitive-state injection trace at a physical boundary, with explicit errors for unsupported treatments. |

The coefficients and injection materials can be prepared outside the
OpenMP cell sweep. The corresponding cell-local laws use values such
as coverage fractions, boundary normals, cell area and physical time,
not a global loop index. Periodic injection windows are evaluated
over the requested physical interval; they are not tied to a retired
Runge–Kutta stage counter.

### `thermal_laws/` — thermal source physics

| File | Responsibility |
| --- | --- |
| `thermal_laws_2d.f90` | Dispatch thermal models through the common treatment interfaces, including scalar-temperature evaluators for thermal-only Newton groups. |
| `thermal_cooling_2d.f90` | Thermal cooling source and its supported explicit, implicit and analytic evaluations. |
| `thermal_cooling_rate.inc` | Shared cooling-rate formula used by the REAL/COMPLEX and full-state/scalar-temperature paths. |

Thermal-only Newton groups can evaluate cooling from a reconstructed
temperature without rebuilding an entire physical-state object.
Combined groups still use the complete physical state where required.
This optimization must not change the underlying cooling law or lose
complex-step derivatives.

### `pore_pressure/` — pore-pressure evolution

| File | Responsibility |
| --- | --- |
| `pore_pressure_2d.f90` | Dispatch pore-pressure source treatments under their standard interfaces. |
| `pore_pressure_darcy_2d.f90` | Darcy-based pore-pressure physics and its supported temporal treatment. |

Darcy is currently the implemented pore-pressure model. Its supported
treatment must be selected through the catalog; an unavailable
interface is not equivalent to a zero physical contribution.

## Architectural boundaries and extension guidance

When adding or modifying physical models:

1. **Keep cell evaluations local.** Provide physical values, source
   parameters and any needed geometric descriptors as arguments. Do
   not pass computational cell indices or traverse `q(j,k,:)` inside
   `physics`. Whole-domain loops and OpenMP scheduling belong outside
   the physical law.
2. **Use canonical state APIs.** Obtain density, temperature,
   composition and other shared properties from `state` rather than
   duplicating their definitions in separate laws.
3. **Preserve all four treatment interfaces.** Implement the
   treatments genuinely supported by the model. Unsupported
   treatments must fail explicitly; do not create compatibility
   wrappers that only forward arguments or report an invented zero
   rate.
4. **Preserve REAL/COMPLEX consistency.** Reuse an `.inc` only for a
   substantive formula that must be the same in both precision paths.
   Do not move elementary local algebra into an include merely to
   remove a short repetition.
5. **Separate evaluation from mutation.** The laws may return
   source rates or candidate substrate changes, but scheduling,
   transactional commit/rollback, and the full time-step acceptance
   decision remain in the numerical layer.
6. **Register new sources explicitly.** Extend
   `numerics/source_catalog_2d.f90`, the corresponding dispatcher
   and configuration/registry mapping as needed. Declare model
   capabilities accurately; do not silently expose unsupported
   Newton or analytic treatments.
7. **Keep physical scratch local or worker-owned.** Prefer
   preallocated reusable workspaces and precomputed immutable
   material descriptors outside performance-critical per-cell
   loops. Avoid hidden mutable shared state in hot evaluation paths.

These boundaries also make the source interfaces suitable places for
future alternative evaluators, including data-driven or ML-assisted
physical closures, provided their contracts, conservation properties,
and numerical treatment requirements are enforced and validated.

## Coding conventions, tests and documentation

All Fortran modules and shared mathematical includes follow
[the mandatory style guide](../../docs/FORTRAN_STYLE.md): four-space
lexical indentation, a maximum of 120 columns for code, trailing
continuation markers aligned to column 120, Doxygen prose within 100
columns, and documentation of procedure arguments and return values.
Main dispatcher/entry-point routines precede auxiliary helpers;
substantial algorithmic blocks should explain invariants and
non-obvious ordering rather than narrating assignments.

Useful checks from the repository root:

```bash
python3 tools/docs/check_fortran_style.py src/
```

```bash
python3 tools/docs/check_fortran_style.py --structural-report src/
```

```bash
python3 tools/docs/check_fortran_style.py --line-width-report src/
```

```bash
python3 TESTS/TEST_SOURCE_DOCUMENTATION/check_documentation.py
```

```bash
python3 .github/scripts/check_fortran_dependencies.py
```

After static checks, compile the solver and execute the relevant
source, state and numerical regression tests, plus representative
examples. Changes to floating-point formulas, branch conditions,
mass/energy accounting or complex-step behavior require scientific
regression evidence, not formatting-only comparison.

The broader architecture and implementation history are described in
[the state architecture roadmap](../../docs/STATE_ARCHITECTURE_REDESIGN_ROADMAP.md).
