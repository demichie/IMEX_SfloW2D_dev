# IMEX-SfloW2D — Development roadmap

**Current development branch:** `feature/two-layer`. This file contains *current
priorities*, not a cumulative commit diary. For authoritative behavior always
inspect the branch HEAD, source code and the CI outcome for that exact revision.

## Current objective: consolidate integrated one- and two-layer execution

Both `N_LAYERS=1` and `N_LAYERS=2` select the standard solver path
without an environment activation gate. The single-layer path remains
the regression baseline. Optimized, strict/debug, OpenMP and Etna
two-layer end-to-end integration all passed on `7388eeeb7`:
[workflow 45822522249](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/45822522249).
This verifies the functional removal of the previous activation gate;
it does not establish general physical validity. The characteristic
restrictions are in [two-layer hyperbolicity](TWO_LAYER_HYPERBOLICITY.md).

The solver-integration workflow now includes the short Etna two-layer
test in every manually dispatched run; extended Etna checks remain
optional. Every changed SHA requires its own CI verification.

The [two-layer architecture and roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
is the single authoritative specification of layer ownership, contacts,
conservative exchanges, spatial coupling, Lie scheduling and scientific
assumptions. Do not introduce competing implementation plans.

| Area | Integrated implementation | Evidence and remaining scope |
| --- | --- | --- |
| Layer-aware state and spatial operator (TL00–TL06) | Canonical `q/qp(component,layer,x,y)`, per-layer reconstructions and boundaries, shared two-layer spatial workset | Structural regressions and optimized/strict-debug integration verified; general coupled hyperbolicity is not assumed. |
| Coupled Lie and physical sources (TL07–TL11) | Local external sources, conservative interlayer mass transfer, and analytic interlayer drag share one scheduler | Source, conservation and wet/dry regressions verified; physical closure calibration remains model-dependent. |
| Input, output and restart (TL12) | Direct `N_LAYERS=1/2` selection, per-layer input and output, shared-substrate budgets and binary continuation | End-to-end Etna checks include output, restart and OpenMP; no special activation is required. |
| Functional integration (TL13a) | **Complete**; short two-layer Etna smoke is now standard solver-workflow coverage | [Workflow 45822522249](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/45822522249) passed all four jobs on `7388eeeb7`. |
| Scientific validity (TL13b) | Formulation retained with documented applicability limits | Inverted density, degenerate equal-density states, excessive shear, coupled spectral behavior and empirical closure calibration are not generally validated. |

## Consolidation scope

The first repository cleanup removes tracked interpreter caches, macOS
metadata and three N7 stage-observer directories whose executable
entrypoints had already been retired. Their scientific contracts remain
named in `TESTS/ACCEPTANCE/acceptance_plan.json` with the status
`historical_only_pending_current_revalidation`; they are recoverable
from [the pre-consolidation revision](https://github.com/demichie/IMEX_SfloW2D_dev/tree/6650582fea44ca5bcbeeb73fc5f631ba0e7bd166/TESTS).
This is **not** new evidence that the missing N7 criteria pass.

Structural rheology consolidation retires the registered Newton-only
demonstration model `rheology_newmodel` and its shared include file,
updates all source IDs/masks and build dependencies, and moves the
independent four-case linear Newton/complex-step check to a fixed-temperature,
zero-yield subcase of the **existing** lava viscous law. The seven
production rheologies' shared REAL/COMPLEX regressions are retained.
All per-law treatment entry points remain in place, including intended
zero contributions and explicit unsupported-method diagnostics.
This removes a development fixture without registering a replacement law.
The change requires optimized, strict/debug, source-physics and integration
Actions on its own commit.

Physical laws, source ordering and numerical results must remain unchanged.
The existing unsuffixed single-layer NAMELISTs are intentional input
contracts and are **not** cleanup targets.

## Invariants of the integrated solver

- **One authoritative state:** `q/qp(component,layer,x,y)`; no synchronized
  copies, global duplicates or transitional internal APIs.
- Keep the single-layer `<FAMILY>_PARAMETERS` input variants alongside the
  explicit `<FAMILY>_L1_PARAMETERS` names where defined. These are intentional
  single-layer input names, **not legacy code**.
- Physics kernels are cell-local. The spatial operator owns pressure/flux
  coupling; the Lie executor owns ordered sources and their limiters.
- An invalid timestep is fatal; no recoverable global substrate rollback,
  retry, post-update clipping or unapproved changes of physics.
- Tests must cover REAL/COMPLEX consistency, conservation, dry states,
  strict/debug and optimized builds, and OpenMP.
- `N_LAYERS=2` is directly runnable, but individual simulations must
  be interpreted within the documented density/shear limits.
  Executable availability is not general scientific validation.

## Current companion references

- [Two-layer hyperbolicity](TWO_LAYER_HYPERBOLICITY.md): characteristic
  limits and accepted scientific restrictions.
- [Two-layer closure audit](TWO_LAYER_CODE_AUDIT.md): historical
  implementation evidence; dated gate statements do not describe the
  current supported execution path.
- [State architecture](STATE_ARCHITECTURE_REDESIGN_ROADMAP.md): current state
  ownership and performance-validation contract.
- [Lie timestep transaction](LIE_TIMESTEP_TRANSACTION_ROADMAP.md): authoritative
  trial/publication and fatal-error semantics.
- [Active-domain OpenMP roadmap](ACTIVE_DOMAIN_OPENMP_ROADMAP.md):
  implemented shared wet/dry update, ping-pong OpenMP halo, fused counting
  and parallel ordered compaction; integration Actions and performance validation pending.
- [NetCDF outputs](NETCDF_OUTPUT.md): visualization schema; binary checkpoint
  is the sole restart format.
- [Fortran style](FORTRAN_STYLE.md): mandatory conventions.
- [Repository documentation index](README.md): guidance and dated reports.

## Archived development history

The former extensive master/Lie/RK migration diary and intermediate plans are
recoverable from the immutable [pre-cleanup repository snapshot](https://github.com/demichie/IMEX_SfloW2D_dev/tree/d91843a0faf8c8491cf50c240cdd50a7cd006407/docs),
including the retired `docs/archive/roadmaps/` directory.
They are **historical evidence, not specifications to reimplement**.
Acceptance report summaries remain under [acceptance](acceptance/README.md),
with links to the raw evidence at the same revision.

**Verification rule:** A green workflow only verifies the SHA it ran.
Structural tests or single-layer simulations alone do not establish
coupled-model validity; functional and scientific evidence remain distinct.
