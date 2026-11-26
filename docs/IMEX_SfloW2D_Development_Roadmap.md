# IMEX-SfloW2D — Current development roadmap

**Baseline:** one- and two-layer execution is integrated on `master`.
This document lists current contracts and outstanding work, **not the
chronology of past migrations**. The source and successful CI on each
specific revision are the authority for verification claims.

## Integrated capabilities

The standard executable accepts `N_LAYERS=1` and `N_LAYERS=2` directly in
`IMEX_SfloW2D.inp` without a special activation environment variable.
The single-layer solver remains the regression baseline.

| Area | Integrated behavior | Limitation |
| --- | --- | --- |
| State | Canonical `q/qp(component,layer,x,y)`, shared terrain and substrate | Preserve dry states, budgets and REAL/COMPLEX agreement |
| Spatial method | Shared active domain, high-resolution fluxes and coupled pressure terms | A CFL step is not proof of two-layer hyperbolicity |
| Source terms | Ordered Lie groups, per-layer external laws, conservative interlayer exchange and analytic drag | Constitutive calibration depends on the scenario |
| Input and output | Per-layer initialization, physical NetCDF fields and one complete binary checkpoint | Current development schemas only |
| Functional integration | Two-layer path in ordinary solver execution and standard short Etna integration test | Engineering integration does not establish general scientific validity |

Earlier optimized, strict/debug, OpenMP and Etna integration results are
recorded by [workflow 45822522249](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/45822522249)
on `7388eeeb7`. They verify **that SHA**, not later commits.

The [two-layer architecture](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md),
[characteristic analysis](TWO_LAYER_HYPERBOLICITY.md) and
[interlayer closure checklist](TL13B_C_PHYSICAL_CLOSURES.md) are the
current detailed scientific references.

## Permanent architecture and numerical contracts

- **One authoritative state:** no duplicated permanent layer arrays,
  synchronized global coefficients, parallel production code paths or
  transitional wrappers.
- **Input names:** unsuffixed single-layer `<FAMILY>_PARAMETERS` NAMELISTs
  are intentional, not compatibility debris.
- **Clear ownership:** `src/state/` owns state and conversions;
  `src/physics/` owns cell-local laws; `src/spatial/` owns fluxes and
  pressure; `src/numerics/` owns source scheduling and accepted updates.
- **Fatal invalid steps:** no recoverable whole-step substrate rollback,
  timestep retry or post-update clipping. Geometry-local deferred
  publication is a distinct mechanism.
- **Necessary temporary data:** trial states, read-only per-layer
  configuration snapshots, pre-CFL primitive information and thread-local
  scratch are required numerical workspaces, not alternate authoritative
  representations. The validated source registry is setup-only; the
  executor owns the operational ordered groups.
- **Scientific integrity:** do not change physics, evaluation order,
  dry-state conventions, source splitting, tolerances or REAL/COMPLEX
  formulas as a by-product of maintenance.

For actual implementation details see the
[state guide](../src/state/README.md),
[active-domain guide](../src/geometry/README.md),
[Lie timestep contract](../src/numerics/README.md) and
[physics-source guide](../src/physics/README.md).

## Outstanding priorities

1. **Two-layer scientific acceptance (TL13b):** examine density/shear
   states and characteristic roots in representative runs; advance
   quantitative closure calibration and the C1–C3/C5 checks in
   [the closure checklist](TL13B_C_PHYSICAL_CLOSURES.md). No universal
   hyperbolicity or empirical validation is implied.
2. **Historical N7 revalidation:** criteria still marked
   `historical_only_pending_current_revalidation` in the
   [acceptance inventory](../TESTS/ACCEPTANCE/acceptance_plan.json)
   need checks using current production observables. Do not restore
   retired stage-observer infrastructure.
3. **Lateral Dirichlet/CFL policy:** independently audit the existing
   production override `boundary_source_enabled=.FALSE.`, rather than
   silently changing it during documentation or structural cleanup.
4. **Local-domain and active-tile architecture (H01–H02):** prepare
   MPI without requiring MPI initially. Introduce explicit local/owned/halo
   extents, then active rectangular tiles, validating numerical equivalence
   to the existing global compact worklists.
5. **MPI and hybrid parallelism (H03–H04):** distribute the owned domain,
   exchange halos, reduce global diagnostics and support MPI+OpenMP with
   distributed output/restart only after H01–H02 acceptance.
6. **OpenMP performance and HPC tuning (H05):** validate full-solver
   reproducibility and measure scaling on matched cases before claiming
   speedups or introducing accelerator-oriented changes.
7. **Distribution and documentation:** publish README, MkDocs, Doxygen
   and GHCR releases only against the deliberately selected, verified
   `master` revision. Publication is manual during development or
   automatic after validation of a stable GitHub Release at `master` HEAD.

## Phase H — Local domains, active tiles and MPI/HPC

This sequence restores the explicit H01–H05 plan from the
[pre-Lie development roadmap](https://github.com/demichie/IMEX_SfloW2D_dev/blob/d91843a0faf8c8491cf50c240cdd50a7cd006407/docs/archive/roadmaps/IMEX_SfloW2D_Development_Roadmap_pre_Lie_2026-10-25.md).
The first two milestones are **serial/OpenMP architectural preparation**;
they must be accepted before MPI domain distribution changes production
numerics. Work on scientific TL13b/N7 acceptance can proceed separately.

### H01 — Local domain and explicit halos

- Define global, locally owned and storage extents, index conversions and
  ownership for cell-centered variables, faces and nodal `B_vertex`.
- Start with `local == global`, retaining the existing scientific
  operators and compact active-cell representation. Make the halo width
  an explicit consequence of actual stencils, **including the two-cell
  reach of the five-point least-squares terrain derivatives**.
- Separate physical boundary filling from future inter-rank exchange.
  Neither physical kernels nor the state representation should depend
  on MPI rank IDs.
- **Acceptance:** equivalence to the current serial solver, including
  HP-PCCU, wet/dry, physical boundaries, LS slopes/Hessian, nodal
  geometry and time-step/conservation checks.

### H02 — Rectangular active tiles

- Partition each **owned local domain** into rectangular tiles, including
  tiles with no active cells; tiles are scheduling and locality units,
  **not independent physical states or alternative solvers**.
- Reuse the authoritative shared wetness mask and workset semantics.
  Preserve canonical cell/face/HP/reconstruction membership and
  deterministic traversal where numerical evaluation order requires it.
  The existing compact worklist is the reference for differential tests.
- Define tile-boundary stencil access and safe ownership of face
  calculations, topographic vertex updates and source/geometry work.
  Do not duplicate global state or add competing persistent cell lists.
- **Acceptance:** multiple tile sizes reproduce the compact-worklist
  reference, including cross-tile neighbors, dry fronts, coupled L1/L2,
  evolving topography and OpenMP 1/2/4/8-thread behavior. Measure locality
  and total-run performance without presupposing a speedup.

### H03 — MPI domain decomposition

- Introduce a 2-D Cartesian rank decomposition over the accepted
  local-domain/active-tile infrastructure. Exchange state and geometric
  halos without physical-variable-specific communication kernels.
- Define **unique or consistently synchronized ownership** of the shared
  nodal `B_vertex`. Assemble contributions from all cells touching a
  vertex, including cells on different ranks, without rank-local
  renormalization. Preserve HP, LS derivative and curvature stencils.
- Use collective reductions only for genuinely global quantities:
  CFL bounds, diagnostics, conservation and termination/failure status.
  The cell-local constitutive laws and ordered Lie groups remain free
  of MPI calls.
- **Acceptance:** MPI 1-rank matches serial; 2/4-rank Cartesian
  decompositions, including nonsquare partitions, reproduce fields,
  boundary flux budgets and restart within stated floating-point limits.

### H04 — Hybrid MPI+OpenMP and parallel I/O

- Validate rank/thread ownership, synchronized halos and OpenMP
  reproducibility. Parallel NetCDF writes only owned interiors and
  uniquely owned vertices; halos are never exported as independent state.
- Define and verify distributed restart for all authoritative fields,
  including substrate/terrain and stochastic state.
- **Acceptance:** supported MPI+OpenMP configurations reproduce reference
  scientific fields and global budgets; distributed output/restart
  reconstructs the same scientific state as serial execution.

### H05 — Performance tuning and accelerator readiness

Only after correctness is verified, profile full-solver scaling and
communication/computation costs, tune tile scheduling and consider
accelerator-friendly kernels. Do not introduce GPU, distributed sparse
backends or new numerical methods as implicit requirements of H01–H04.

**Unchanged constraints:** one authoritative state and physical model;
no duplicate tile-local copies of complete fields; source ordering,
REAL/COMPLEX behavior, positivity/admissibility, conservation, source
limits and wet/dry behavior must remain correct. Any change to physics
or stencil width requires a separately justified scientific decision.

## Verification and history

Use the three independently dispatchable workflows for source/documentation,
physics/numerical regressions and solver/OpenMP integration. Published stable
releases invoke them automatically before distribution; include
optimized and strict/debug checks, one- and two-layer examples, binary
restart and conservation tests as applicable. **Green CI verifies only
the SHA on which it ran.**

- [Workflow catalogue](../.github/ACTIONS.md)
- [Scientific acceptance inventory](../TESTS/ACCEPTANCE/acceptance_plan.json)
- [Dated scientific evidence](acceptance/README.md)
- [Superseded development notes at a fixed Git revision](https://github.com/demichie/IMEX_SfloW2D_dev/tree/cccdee3ef8d54b7efea92ace97d4beabc0c52ac5/docs)

**Functional integration is complete; unrestricted physical approval is not.**
