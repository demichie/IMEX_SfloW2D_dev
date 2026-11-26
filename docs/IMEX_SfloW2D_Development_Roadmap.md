# IMEX-SfloW2D — Development roadmap

**Current development branch:** `feature/two-layer`. This file contains *current
priorities*, not a cumulative commit diary. For authoritative behavior always
inspect the branch HEAD, source code and the CI outcome for that exact revision.

## Current objective: TL13

The single-layer path remains the regression baseline. The integrated
two-layer architecture can now run directly with `N_LAYERS=2`; no
environment activation gate remains. Functional integration passed the
previously recorded CI checks; scientific validity is restricted by the
conditions in [two-layer hyperbolicity](TWO_LAYER_HYPERBOLICITY.md).
Direct execution is not a claim of general physical validation.

The [two-layer architecture and roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
is the single authoritative specification of layer ownership, contacts,
conservative exchanges, spatial coupling, Lie scheduling and scientific
assumptions. Do not introduce competing implementation plans.

| Stage | Current position | Next required evidence |
| --- | --- | --- |
| TL00–TL06 | Shared layer-aware state, spatial reconstruction/BC and coupled spatial fixtures implemented | Preserve single-layer output; confirm shared wet/dry workset and per-layer boundaries in end-to-end cases |
| TL07 | Coupled Lie trial, source scheduling and optional layer physics exercised structurally | Validate actual two-layer production timestep |
| TL08–TL11 | Interlayer settling, elutriation, carrier and analytic drag implemented/tested in isolated and integrated fixtures | Validate coupled mass/species/energy/momentum budgets through a complete run |
| TL12 | Two-layer input, snapshots, named NetCDF physical fields, and binary restart structurally implemented | Verify real driver I/O, restart continuation and OpenMP |
| TL13a | **Complete — functional integration verified on `89a9f2b`** | Run `44712747883` verified optimized and strict/debug integration, OpenMP and two-layer Etna end-to-end tests; the no-gate integration change needs new CI evidence. |
| TL13b | **Open** | Scientifically validate hyperbolicity/shear, density contrasts, coupled pressure, wet/dry, convergence and physical closures |

## Invariants preserved at TL13a closure

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
- [Two-layer closure audit](TWO_LAYER_CODE_AUDIT.md): technical baseline and
  remaining functional evidence; historical commit-specific claims must be
  rechecked against the current HEAD.
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
