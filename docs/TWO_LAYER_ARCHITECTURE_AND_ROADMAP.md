# Two-layer architecture and implementation status

The solver supports `N_LAYERS=1` and `N_LAYERS=2` directly from
`IMEX_SfloW2D.inp`, without any environment activation flag. The
layer-aware architecture is integrated in the standard solver, not
an alternate development implementation. The existing single-layer
numerical path remains the regression baseline.

**Scientific scope:** the coupled model has explicitly accepted
applicability limits. Consult
[TWO_LAYER_HYPERBOLICITY.md](TWO_LAYER_HYPERBOLICITY.md) for the
characteristic spectrum and CFL limitations, and
[TL13B_C_PHYSICAL_CLOSURES.md](TL13B_C_PHYSICAL_CLOSURES.md)
for mass-transfer, drag and energy-budget assumptions. A completed
run is not a certificate of hyperbolicity or empirical calibration.

This is the **current architecture specification**, not a chronological
implementation report. The dated TL13a development audit is retained at
[a fixed historical Git revision](https://github.com/demichie/IMEX_SfloW2D_dev/blob/cccdee3ef8d54b7efea92ace97d4beabc0c52ac5/docs/TWO_LAYER_CODE_AUDIT.md);
it does not describe the present runtime activation status.

## 1. Canonical ownership and supported configurations

- **State:** one canonical `q/qp(component,layer,x,y)` representation.
  Model layout owns the equation and primitive-component indices;
  no duplicated layer globals, synchronized copies or transitional
  state interfaces are maintained.
- **Per-layer configuration:** rheology, thermal behavior,
  pore pressure/Darcy, stochastic/OU and external mass-exchange
  sources. The species properties, geometry, atmosphere, terrain,
  deposition and erodible-substrate inventories are shared.
- **Contacts:** `layer_contacts_2d` classifies wet/wet, L1-only,
  L2-only and dry/dry cells from the current candidate. The
  uppermost wet layer contacts atmosphere, and the lowermost wet
  layer contacts terrain. These contacts are derived, not copied
  into persistent duplicated configuration.
- **Cold starts:** settling can populate an initially dry L1;
  elutriation can populate an initially dry L2. Interlayer drag
  requires both layers physically present.
- **Input:** the single-layer forms `<FAMILY>_PARAMETERS` and
  `<FAMILY>_L1_PARAMETERS` have intentional single-layer semantics,
  not legacy compatibility. The two-layer case requires explicit
  L1/L2 input groups. `INTERLAYER_PARAMETERS` is the unique
  pair-wide configuration for mass and momentum transfer.
- **Domain and boundaries:** both layers use the same terrain,
  geometry and shared OR active workset, with per-layer boundary
  data on all four sides. Inactive state components remain zero.

`src/physics/` owns cell-local physical laws; `src/spatial/`
owns reconstruction, numerical fluxes and nonconservative pressure;
`src/numerics/` owns Lie ordering, candidate checks and publication;
`src/state/` owns state and conversions. Physical kernels do not
own grid-cell indices or duplicate the timestep integrator.

## 2. Coupled spatial pressure and characteristic bounds

The layer surfaces and gravity metric are derived from the
canonical terrain `B`:

```text
eta1 = B+h1
eta2 = B+h1+h2
G = 1/(1+B_x**2+B_y**2)
```

With `L=Gamma2*h2=rho2*gprime2*h2`, the implemented interlayer
pressure contributions in horizontal direction `j` are

```text
S_inter1 = -G*h1*d_j(L)
S_inter2 = -L*d_j(G*h1).
```

The remaining L2 pressure contribution is already included in
the hydrostatic path and must not be duplicated. The existing
nonnegative L2 buoyancy-pressure clamp is retained. Interlayer
pressure is spatial, **not** a Lie source. The geometry metric
and curvature derive from `B`, not `B+h1`.

The coupled spatial operator reuses both layer reconstructions
and the shared face workset. Face speeds use the existing
**total-depth envelope**: each present layer contributes its
single-layer bounds evaluated at `H=h1+h2`. For smooth,
positive-buoyancy, stably stratified wet/wet states this
encloses the **real** gravitational characteristic speeds,
under the assumptions in
[TWO_LAYER_HYPERBOLICITY.md](TWO_LAYER_HYPERBOLICITY.md).
It is not a coupled eigensystem or a hyperbolicity detector;
complex roots cannot be cured by choosing a smaller timestep.

## 3. Lie integration and external source ownership

The solver has one coupled timestep and one cell-outer/OpenMP
Lie executor. Groups are ordered by numeric ID; within the same
ID the scopes are ordered `L1`, `L2`, then `INTERLAYER`.
The treatment method of each source is constrained by the
canonical source registry.

**Group 0** is evaluated on the pre-spatial state, but its
limited increment is applied to the post-spatial candidate:

```text
q*       = q^n + dt*L(q^n)
R0       = S0(q^n, qp^n, ...)
R0_limit = Limit(R0, q*, dt, inventories)
q^(0)    = q* + dt*R0_limit
```

Subsequent groups read accepted local candidates. Interlayer
transactions evaluate both layers from one shared snapshot,
limit simultaneous donor withdrawals, and commit paired sources.
Admissibility is checked before publication. There is no
after-the-fact mass clipping, recoverable global rollback or
retry of rejected timesteps; invalid candidates terminate with
diagnostics.

Layer-local mass exchange with substrate and atmosphere is
independent of L1/L2 transfer. **Static substrate deposition
is restricted to L1**. L2 particles first settle into L1;
an isolated wet L1 may entrain ambient air, whereas an isolated
wet L2 may contact terrain where a configured physical law permits.

The two directions of interlayer mass transfer share one
**EXPLICIT** Lie group. Chézy interlayer drag is separately
integrated in a positive-ID **ANALYTIC** group. Source group
ordering is configurable and affects the evolving state,
but it does not change the spatial characteristic polynomial.

## 4. Physical exchange and conservation

`src/physics/interlayer_mass_exchange/` implements independent
proposed per-species volume rates: settling from L2 to L1 and
elutriation from L1 to L2. The coordinator constrains the
**gross** solid and carrier withdrawals against each donor's
original substep inventory, both final packing states and
conservative admissibility. Separate limit factors can apply
to settling, elutriation and, where allowed, associated carrier
flow. Newly received material is not automatically spent in
the same Lie transaction.

`ELUTRIATION_SOLID_VOLUME_FRACTION` is the fraction **in the
transferred stream**, not the target composition of L2:

- `-1`: derive it from the current L1 donor inventory;
- a prescribed `0<alpha_lift<=1`: use the specified value.

For accepted lifted solid volume `dV_s`, requested donor
carrier volume is

```text
dV_c = (1-alpha_lift)/alpha_lift*dV_s.
```

With `ELUTRIATION_CARRIER_MODE='LIMIT_SOLIDS'`, insufficient
carrier reduces both transfers proportionally. With
`'ALLOW_CONCENTRATED'`, less carrier may accompany the
accepted solid flow, producing a more concentrated lifted
stream without creating gas or violating final packing.
Atmospheric entrainment into L2 is an **external** process
that can later dilute the upper mixture.

Paired transfers conserve species/carrier **mass** between
layers. Carrier **volume** may differ between donor and
recipient when gas densities differ. Momentum and sensible
thermal transfers use the applicable donor/interfacial
velocities and thermodynamic properties. The physical-closure
audit documents the exact budgets and empirical uncertainties.

Chézy drag conserves combined momentum at fixed layer
masses and reduces relative kinetic energy. It does not
return the dissipated mechanical energy to sensible heat.
**This no-heating approximation is accepted**, not an
unimplemented feature.

## 5. Input, output, restart and examples

The same executable and input parser handle one or two layers.
The two-layer case uses explicit layer initial conditions,
per-layer snapshots and named NetCDF physical fields, a
common terrain/substrate state, and layer/total budgets.
A single canonical binary checkpoint contains the
selected layer count and full restart state; no obsolete
versioned checkpoint readers are retained.

Executable examples and their scientific interpretations:

- [Etna granular collapse and dilute upper current](../EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/README.md):
  L1 elutriation, L2 air entrainment, external L1 deposition
  and optional interlayer drag.
- [One-dimensional two-layer dam-break](../EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D/README.md):
  spatial transport with initially dry L2 and no transfers.

See [NETCDF_OUTPUT.md](NETCDF_OUTPUT.md) and
[RESTART_INPUT.md](RESTART_INPUT.md) for the file contracts.

## 6. Verification and model acceptance

The baseline integration run
[45822522249](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/45822522249)
passed optimized, strict/debug, end-to-end OpenMP comparison,
and Etna two-layer integration on commit `7388eeeb7`.
Focused state, physical-source and spatial regressions were
also verified on earlier documented code revisions.

`solver-validation.yml` now runs the short Etna two-layer
integration as a **mandatory** job on manual and release dispatch.
Extended Etna conservation/restart/OpenMP tests remain optional
manually but are enabled automatically for a stable published
release. Both standard integration profiles remain enabled.
Ordinary pushes do not trigger Actions; every change must be
validated on its **own SHA**.

**Functional acceptance:** the one- and two-layer execution
paths are integrated, with no special activation variable
and no redundant two-layer solver architecture.

**Scientific limitations:** density inversion has complex
gravity-wave roots already at equal velocities, equal
densities give a degenerate gravitational block, and
sufficient interlayer shear may cause complex roots even
with stable density ordering. The total-depth CFL envelope
does not establish hyperbolicity. Dry states, nonpositive
buoyancy, variable composition and empirical transfer
closures carry additional documented restrictions.
No general coupled-model physical calibration or strong
hyperbolicity for all admissible states is claimed.

No new physical law, pressure correction, hyperbolicity
limiter, timestep workaround or new activation gate is
introduced by this engineering completion. The historical
TL13a implementation tests and TL13b scientific decisions
remain documented by the linked audits, not as runtime
features to reactivate.
