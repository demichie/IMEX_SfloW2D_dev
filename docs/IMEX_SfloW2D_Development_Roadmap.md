# IMEX-SfloW2D — Development roadmap

**Current development branch:** `feature/two-layer`. This file contains *current
priorities*, not a cumulative commit diary. For authoritative behavior always
inspect the branch HEAD, source code and the CI outcome for that exact revision.

## Current objective: TL13

The single-layer production path remains the reference baseline. The two-layer
architecture is structurally implemented but **`N_LAYERS=2` remains disabled
in the production driver** until the controlled end-to-end functional gate and
subsequent scientific validation.

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
| TL13a | **Open** | Controlled two-layer end-to-end simulation with full parsing, IC/BC, evolving L1/L2, budgets, output, restart and determinism |
| TL13b | **Open** | Scientifically validate hyperbolicity/shear, density contrasts, coupled pressure, wet/dry, convergence and physical closures |

## Invariants not to change during TL13a

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
- Do not enable `N_LAYERS=2` as a validated production configuration before
  the functional and scientific gates are independently satisfied.

## Current companion references

- [Two-layer closure audit](TWO_LAYER_CODE_AUDIT.md): technical baseline and
  remaining functional evidence; historical commit-specific claims must be
  rechecked against the current HEAD.
- [State architecture](STATE_ARCHITECTURE_REDESIGN_ROADMAP.md): current state
  ownership and performance-validation contract.
- [Lie timestep transaction](LIE_TIMESTEP_TRANSACTION_ROADMAP.md): authoritative
  trial/publication and fatal-error semantics.
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
Structural two-layer tests or single-layer production simulations do not imply
TL13a completion, and TL13a does not imply TL13b scientific approval.
