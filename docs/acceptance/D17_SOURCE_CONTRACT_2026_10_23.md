# D17 source-term inventory and compatibility contract — 23 October 2026

Status: **inventory and contracts documented; runtime D17 not implemented**.
Baseline: `1310511ed940dc50325f656236ca7ce467eb8d52` (includes the recorded D16 Intel ifx CI acceptance).

Historical phase-A specification. The subsequent [D17-B implementation and
validation report](D17_B_SOURCE_REGISTRY_2026_10_23.md) records the runtime metadata
now present; the phase-A contracts and baseline measurements below remain intact.
No production source, input, workflow, tolerance or numerical algorithm changes
are part of this phase. No simulation or new compiler/CI acceptance is claimed.

The [machine-readable manifest](D17_SOURCE_TERMS_MANIFEST_2026_10_23.json)
contains 34 entries: 30 physical contributions and four dependent equation
assemblies. `SOURCE_TERM_*` names below are proposed identifiers, not existing
Fortran enums. They inventory all nonzero rheology rate branches, cooling,
pressure diffusion, cell injections, boundary-state injection and mass exchange.
Activation is conditional: mutually exclusive rheologies are not all registered
in one simulation. Hashes and source anchors identify the reviewed baseline.

## Default groups and scheduling

Group identity is `(cadence,id)`. Method is homogeneous within that identity;
ID alone does not determine execution phase or order. A proposed explicit stage
group **4** avoids colliding with existing implicit IDs 1–3. Stage group 1 and
step group 1 remain distinct. Inactive energy/pressure groups may be omitted;
internal dense indexing must not overwrite user-facing IDs or execution order.

| Cadence / ID | Method | Target | Existing behavior to preserve |
|---|---|---|---|
| stage / 1 | analytic | temperature / Q4 | Lava cooling first when active |
| stage / 2 | analytic | excess pressure / Qp | Pressure diffusion after the updated temperature |
| stage / 3 | newton | velocity / Q2,Q3 | SI correction then implicit drag and dissipative correction |
| stage / 4 | explicit | cell rates and boundary-state hook | Existing requested EX evaluation and spatial boundary reconstruction |
| step / 1 | explicit | mass-exchange bundle | One post-RK proposal/application, inventories and bed update |

The momentum partition remains present even with rheology disabled: its zero
rate and numerical projections must not disappear merely because no friction
term is active. The source registry does not redefine the spatial explicit
partition or move hydrostatic/bed-slope path terms into cell sources.

The actual time-step order is:

1. Choose the final `dt`; prepare exact whole-step injection-window averages
   with `set_source_time_coefficients(t,dt)`; perform existing OU preparation.
2. For each RK stage, construct the known state from previously stored EX/IM
   and spatial contributions; retain current clipping and state conversion.
3. Apply ordered implicit groups when the diagonal is nonzero and the cell is
   wet. Defaults are temperature, pressure, velocity. In a Newton momentum
   group, evaluate SI on the entry state, correct stopping **before** Newton,
   solve its implicit block, then apply the existing dissipative correction.
4. Form the effective IM contribution from the corrected stage increment.
   Evaluate requested EX terms on the solved stage and evaluate the spatial
   operator with its existing boundary-state hook. Preserve `need_explicit_stage`,
   the tableau weights and the existing `t` argument; do not introduce stage-time
   shifts or extra source evaluations as part of the registry refactor.
5. Complete the RK weighted update, its positivity/dry-state projections and
   primitive reconstruction. These checks already occur before mass exchange;
   a generic new end-of-step checklist must not postpone them.
6. Evaluate the post-RK mass-exchange proposals on one common cell state. Apply
   the source-coverage mask to equation rates, erosion/deposition and bed rate;
   update conservative state and inventories once, perform existing cell checks
   and primitive reconstruction, then project/apply accumulated bed proposals.
7. Advance runtime time and retain existing diagnostics/output scheduling.

## Term inventory

`SI` means a real pre-Newton rate, not a fourth configurable method. `I` means
a real/complex implicit rate. `A/I` has an analytic default with an existing
implicit-rate path for configured Newton groups. `B` is a primitive boundary
state consumed by the spatial flux, not a cell rate. `D` is dependent assembly
and is never independently scheduled. `M` is the common post-RK mass batch.
Equation indices follow the configured layout; optional components are enabled
only when present. Exact activation, state provenance, equation footprints,
current interfaces, owners and evaluation batches are in the manifest.

| Proposed term (without `SOURCE_TERM_`) | Current kernel | Role / default group |
|---|---|---|
| `RHEOLOGY_VOELLMY_COULOMB` | [eval_voellmy_semi_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_voellmy_2d.f90#L189) | SI / stage 3 |
| `RHEOLOGY_VOELLMY_TURBULENT_DRAG` | [eval_voellmy_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_voellmy_2d.f90#L135) | I / stage 3 |
| `RHEOLOGY_PLASTIC_YIELD` | [eval_plastic_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_plastic_2d.f90#L81) | I / stage 3 |
| `RHEOLOGY_LAVA_VISCOUS` | [eval_lava_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_lava_2d.f90#L104) | I / stage 3 |
| `RHEOLOGY_LAVA_YIELD` | [eval_lava_semi_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_lava_2d.f90#L140) | SI / stage 3 |
| `RHEOLOGY_LAHAR_VISCOUS` | [eval_lahar_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_lahar_2d.f90#L146) | I / stage 3 |
| `RHEOLOGY_LAHAR_TURBULENT_DRAG` | [eval_lahar_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_lahar_2d.f90#L146) | I / stage 3 |
| `RHEOLOGY_LAHAR_YIELD` | [eval_lahar_semi_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_lahar_2d.f90#L204) | SI / stage 3 |
| `RHEOLOGY_EMPIRICAL_DRAG` | [eval_empirical_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_empirical_2d.f90#L77) | I / stage 3 |
| `RHEOLOGY_QUADRATIC_DRAG` | [eval_quadratic_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_quadratic_2d.f90#L79) | I / stage 3 |
| `RHEOLOGY_FROUDE_FRICTION` | [eval_froude_momentum_rate](../../src/physics/rheology_laws/rheology_froude_2d.f90#L112) | SI / stage 3 |
| `RHEOLOGY_LUCAS_FRICTION` | [eval_lucas_momentum_rate](../../src/physics/rheology_laws/rheology_lucas_2d.f90#L109) | SI / stage 3 |
| `RHEOLOGY_MUI_FRICTION` | [eval_mui_momentum_rate](../../src/physics/rheology_laws/rheology_mui_2d.f90#L139) | SI / stage 3 |
| `RHEOLOGY_MUI_VOELLMY_FRICTION` | [eval_mui_voellmy_momentum_rate](../../src/physics/rheology_laws/rheology_mui_voellmy_2d.f90#L209) | SI / stage 3 |
| `RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG` | [eval_mui_voellmy_implicit_momentum_rate](../../src/physics/rheology_laws/rheology_mui_voellmy_2d.f90#L157) | I / stage 3 |
| `THERMAL_RADIATIVE_COOLING` | [eval_rheology_cooling_energy_rate](../../src/physics/thermal_laws_2d.f90#L98) | A/I / stage 1 |
| `THERMAL_CONVECTIVE_COOLING` | [eval_rheology_cooling_energy_rate](../../src/physics/thermal_laws_2d.f90#L98) | A/I / stage 1 |
| `PORE_PRESSURE_DISSIPATION` | [eval_pore_pressure_diffusion_source](../../src/physics/pore_pressure_laws_2d.f90#L56) | A/I / stage 2 |
| `BOTTOM_CIRCULAR_INJECTION` | [eval_expl_terms](../../src/physics/equation_terms_2d.f90#L109) | EX / stage 4 |
| `BOTTOM_FISSURAL_INJECTION` | [eval_expl_terms](../../src/physics/equation_terms_2d.f90#L109) | EX / stage 4 |
| `RADIAL_INJECTION` | [eval_expl_terms](../../src/physics/equation_terms_2d.f90#L109) | EX / stage 4 |
| `GEOMETRY_CURVATURE` | [eval_curvature_momentum_source](../../src/physics/equation_terms_2d.f90#L1517) | EX / stage 4 |
| `BOUNDARY_INJECTION` | [eval_source_bdry](../../src/physics/equation_terms_2d.f90#L1169) | B / stage 4 |
| `MASS_EXCHANGE_EROSION` | [eval_erosion_rates](../../src/physics/mass_exchange_laws_2d.f90#L235) | M / step 1 |
| `MASS_EXCHANGE_DEPOSITION` | [eval_solid_deposition_rate](../../src/physics/mass_exchange_laws_2d.f90#L146) | M / step 1 |
| `MASS_EXCHANGE_CARRIER_LIQUID_LOSS` | [eval_mass_exchange_terms](../../src/physics/equation_terms_2d.f90#L843) | M / step 1 |
| `MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS` | [eval_mass_exchange_terms](../../src/physics/equation_terms_2d.f90#L843) | M / step 1 |
| `MASS_EXCHANGE_CARRIER_DARCY_LOSS` | [eval_pore_pressure_gas_loss](../../src/physics/pore_pressure_laws_2d.f90#L200) | M / step 1 |
| `MASS_EXCHANGE_ENTRAINMENT` | [eval_air_entrainment_rate](../../src/physics/mass_exchange_laws_2d.f90#L304) | M / step 1 |
| `MASS_EXCHANGE_VAPORIZATION` | [eval_vaporization_mass_rate](../../src/physics/thermal_laws_2d.f90#L294) | M / step 1 |
| `INLET_THERMAL_COUPLING` | [eval_inlet_thermal_energy_source](../../src/physics/thermal_laws_2d.f90#L46) | D / parent |
| `INLET_PORE_COUPLING` | [eval_pore_pressure_inlet_source](../../src/physics/pore_pressure_laws_2d.f90#L329) | D / parent |
| `MASS_PORE_COUPLING` | [eval_pore_pressure_mass_exchange_source](../../src/physics/pore_pressure_laws_2d.f90#L355) | D / parent |
| `MASS_STOCHASTIC_COUPLING` | [eval_mass_exchange_terms](../../src/physics/equation_terms_2d.f90#L843) | D / parent |

## Numerical and API contracts

### Analytic results and target blocks

`eval_analytic_terms(group, solver_id, q_before, context, interval, q_after)`
returns a **full conservative state**, initially copied from `q_before`, with
only the declared target block replaced. It returns neither an instantaneous
rate nor an increment. A target is temperature (`Q4`), pore pressure (`Qp`) or
the two-component velocity/momentum block (`Q2,Q3`). A single velocity solve
updates both components together and preserves the current directional/stopping
policy. The assembly wrapper converts primitive results to conservative ones:
`Q4=C_area*T_after`, `Qp=Q1*p_after`, `Q(2:3)=Q1*velocity_after`.

The numerical coordinator passes the effective interval `dt*a_ii` for stage
solves and forms `IM_effective=(q_after-q_known)/(dt*a_ii)` only when the diagonal
is nonzero, after all existing corrections/projections. Zero diagonal stages
follow the existing direct path without division. Group results compose by
replacing the declared target block in the current ordered state, never by
resetting unrelated components to the original step state.

Analytic means the existing frozen-coefficient **exact-flow** law. It does not
mean an algebraic solution of backward Euler and is not interchangeable with
Newton. A physical analytic kernel may receive an explicit interval but knows
nothing about tableau indices, group IDs, cadence or the runtime clock. Ordinary
physical rates remain rates; existing admissibility caps that take an interval
and inventory remain explicit bounded-proposal APIs until separately migrated.

### Newton coupling and dissipative corrections

A solve descriptor owns an ordered list of target equations, associated term
IDs, method and execution order. A Newton block may contain both momentum
components, energy and pressure. All cross derivatives within a legacy joint
block remain in its Jacobian. Splitting a shared block into sequential velocity,
thermal and pressure solves changes the method and is forbidden in compatibility
mode. Multiple genuinely independent blocks may be ordered explicitly; outside
block values are frozen from the current ordered state. Unsupported overlaps or
ambiguous ownership are rejected before integration.

Keep normalized complex-step differentiation, residual bounds, precision-specific
regularization and failure/rollback diagnostics. All selected implicit rates
must reach both residual and Jacobian through the same group/solver filter.
For the default Newton momentum path preserve SI evaluation, pre-Newton stopping,
implicit solve, post-Newton correction and effective-rate assembly. Their
implementation is numerical application policy, not a new physical `term_id`.
The corrected pre-Newton placement from `9b2a4989` is a regression boundary.

### Shared kernels and conserved assembly

A physical ID does not imply a separately callable kernel already exists.
Lahar viscous/turbulent drag share a summed kernel; radiation/convection share
one cooling kernel/analytic solver. Preserve those evaluation and accumulation
batches while introducing metadata. The current analytic velocity solver consumes
all active friction parts of the selected law once. A partial analytic bundle is
unsupported unless the kernel explicitly implements and validates that subset.

Mass exchange is likewise one bounded proposal/application transaction initially:

- Erosion and deposition use the same pre-application state and their existing
  `dt`/inventory caps; deposition supplies porosity loss and vaporization.
- Configured liquid loss, porosity loss and Darcy gas loss are summed **before**
  the shared packing/carrier-reserve limiter. Limiting each term independently
  is not equivalent. Keep the existing summation order and carrier coefficients.
- Entrainment and vaporization enter conservative assembly before dependent
  pressure/stochastic mass couplings. Heat and phase-component contributions are
  dependent parts of these physical bundles, not separately integrated sources.
- Preserve `alphastot_min`, dry-state guards, one source-coverage mask, inventory
  updates, re-erodible deposits, positivity checks and topography projection.

New configurations that split an atomic evaluation batch across methods/cadences
or independently schedule dependent children must be rejected until a separate
conservative contract exists. No physics formula or global floor is introduced
merely to make a proposed grouping expressible.

### Cell rates, boundary states and context

The source catalogue needs a representation/capability descriptor in addition
to method and cadence: cell rate, boundary primitive state, or dependent assembly.
This is not an independent temporal `treatment` selector. Boundary injection
remains `eval_source_bdry -> reconstruction -> numerical spatial flux`; the new
registry selects the hook, never adds its state vector to a conservative rate.
The four proposed generic rate/analytic interfaces do not replace that existing
boundary-state interface. Hydrostatic/bed-slope path forcing remains excluded.

A read-only evaluation context must carry state/geometry, slope and Hessian,
gravity multiplier, frozen local stochastic value, injection coverage/arc
geometry, source-window snapshot and the existing time argument as needed.
Mass-exchange proposals additionally receive interval and inventory snapshots.
Numerical application owns mutations. Define these arguments before coding new
wrappers; do not silently re-read a different global state or draw new noise
inside a term/Jacobian evaluation. Per-step caches have one initializer before
RK and read-only use inside parallel cell loops; per-cell work is local.

Real/complex support is a per-interface capability, not an obligation to
complexify boundary, explicit or analytic paths that are real today. Common
mathematics can use an existing complex kernel with a real wrapper, generic
procedures or narrow `.inc` kernels where useful. Do not create duplicate
formulas or move already extracted D16 laws solely to impose a file pattern.

## Compatibility of legacy input

Continue accepting the current numerical-parameter fields:
`IMEX_VELOCITY_GROUP`, `IMEX_ENERGY_GROUP`, `IMEX_PORE_PRESSURE_GROUP` and
`IMEX_GROUP_METHOD(1:3)`, including case handling and valid reordered/shared
groups. Omitted input means energy=1/analytic, pressure=2/analytic,
velocity=3/Newton. Retain existing namelist backup/restart compatibility and
record the effective normalized mapping for diagnostics.

Normalize legacy configurations once at initialization:

1. Resolve active equation masks and preserve original user IDs before dense
   compaction. Sort implicit work by configured group ID, as today.
2. Keep each shared NEWTON group as **one coupled solve block**, including SI
   only when it contains momentum. Do not replace it with three scalar solves.
3. Expand each shared ANALYTIC group into one internal group per active solve,
   in the actual current order **velocity, temperature, pressure**. Refresh the
   primitive state between solves. This preserves the mixed analytic group that
   is already covered by `test_grouped_implicit_stage.f90`.
4. Allocate fresh internal IDs for expanded groups (reserve stage ID 4 for EX;
   allocate split IDs above 4), retain their original legacy ID and assign an
   explicit execution-order key. Internal ID magnitude is not the execution
   order; existing public input IDs and backup fields remain unchanged.
5. New native source-group input follows the one-analytic-solve rule. Legacy
   mixed analytic input remains accepted through normalization. If both new and
   explicitly supplied legacy source-group configuration are present, reject
   the ambiguity; parser presence tracking must distinguish explicit fields
   from default initialization.

The exact new namelist syntax remains a runtime-phase deliverable; no new input
field or enum is implemented by this document. Validate method/cadence/group,
capability, target and ordering before starting. Register each active physical
contribution once globally across cadences; `stage/step` namespaces do not
license counting the same physical process twice. A zero coefficient is not a
second term. Optional/inactive equations do not authorize out-of-range accesses.

## Implementation gates and verification

| Phase | Deliverable | Status |
|---|---|---|
| A | Baseline catalogue, defaults and compatibility contracts | Documented here; static verification only |
| B | Enum/descriptor modules, immutable registries and legacy normalization | Not implemented |
| C | Analytic assembly bridge and physics-free numerical solver | Not implemented |
| D | Term/solver-filtered real and complex rate interfaces; coupled-block preservation | Not implemented |
| E | EX/boundary and atomic step mass-exchange integration | Not implemented |
| F | New parser/capability validation, regressions and acceptance | Not implemented |

Gate B must test default IDs without method collisions, sparse active groups,
legacy reordering/shared groups, atomic batch membership, duplicate terms across
cadences and invalid target/capability combinations. A source group type needs
an execution-order key, target/block descriptors and representation/capability
metadata in addition to the minimum fields in the roadmap.

Before closing C–F, preserve the existing grouped implicit tests (including
joint Newton residuals, shared analytic energy/pressure, sequential state refresh
and effective stage increments), analytic velocity/cooling/pressure oracles,
SI stopping and the D16 slowdown regression. Add filtered residual/Jacobian
checks, zero diagonal/dry/disabled cases, native-vs-legacy normalized equivalence,
ambiguous input rejection and unsupported partial-bundle rejection. Verify all
injections including curvature and boundary flux without double application;
window crossings, mass conservation, source masking, inventories, positivity,
re-erodible deposits, topography and restart must retain their current behavior.
Use `n_RK=2` for new operational fixtures, `h_x` for lava and `froude_blend` for
other cases. Keep the existing temporal-order tests and their supported tableaux.

Runtime acceptance requires appropriate optimized/strict local tests and the
project's established portability gates. The ifx evidence recorded by baseline
`1310511e` belongs to D16, not to future D17 code. This phase ran only static
checks: manifest syntax, unique IDs, references, group-method consistency,
source hashes, Markdown links and unchanged production/input files. No numerical
test was run because runtime behavior was not changed. D17 remains **IN CORSO**.
