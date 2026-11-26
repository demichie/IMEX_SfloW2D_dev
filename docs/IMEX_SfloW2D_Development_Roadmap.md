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
4. **OpenMP reproducibility and performance:** the shared wetness mask,
   ping-pong halo, parallel ordered worklist and dependent workset
   algorithms are integrated. Validate full-solver behavior and measure
   complete-run scaling on matched cases before claiming speedups.
5. **Distribution and documentation:** publish README, MkDocs, Doxygen
   and GHCR releases only against the deliberately selected, verified
   `master` revision. Publishing workflows require manual dispatch.

## Verification and history

Use the three independent manual workflows for source/documentation,
physics/numerical regressions and solver/OpenMP integration; include
optimized and strict/debug checks, one- and two-layer examples, binary
restart and conservation tests as applicable. **Green CI verifies only
the SHA on which it ran.**

- [Workflow catalogue](../.github/ACTIONS.md)
- [Scientific acceptance inventory](../TESTS/ACCEPTANCE/acceptance_plan.json)
- [Dated scientific evidence](acceptance/README.md)
- [Superseded development notes at a fixed Git revision](https://github.com/demichie/IMEX_SfloW2D_dev/tree/cccdee3ef8d54b7efea92ace97d4beabc0c52ac5/docs)

**Functional integration is complete; unrestricted physical approval is not.**
