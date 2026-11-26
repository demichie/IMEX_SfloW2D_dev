# State architecture — current ownership and verification

**Scope:** implementation contract for state layout, physical conversion and
performance. Superseded migration history is available in the
[pre-cleanup Git snapshot](https://github.com/demichie/IMEX_SfloW2D_dev/blob/d91843a0faf8c8491cf50c240cdd50a7cd006407/docs/STATE_ARCHITECTURE_REDESIGN_ROADMAP.md).
The previous title is retained so existing links from `src/state/README.md`
and `src/physics/README.md` continue to work.

## Canonical source of truth

- `src/state/` contains `state_types`, `state_fields`,
  `state_composition`, `state_mixture`, `state_conversions`,
  `state_physical` and the `state` facade. The former `state_2d`,
  `state_conversion_2d` and `src/state_new/` are not parallel
  production architectures.
- The canonical conservative and primitive storage is
  `q/qp(component,layer,x,y)`. The two-layer layout is an extension of
  these arrays, not an additional collection of L1/L2 state copies.
- Shared topography, geometry, substrate and material inventories are not
  duplicated per layer; layer-specific configurations and optional state
  components belong to the corresponding layer.
- Distinguish conservative equations (mass, momenta, energy, components)
  from primitive depth, velocity, temperature and composition, even when
  their integer indices happen to agree.

## Scientific invariants

1. Preserve dimensional consistency, density/composition closure, dry and
   near-dry conversions, component activation, admissibility and
   conservative budgets. Preserve REAL/COMPLEX mathematical equivalence
   wherever the numerical methods require it.
2. Keep exactly one implementation of each mathematical law, reusing
   existing shared kernels where equivalent. No global coefficient mirrors,
   synchronized temporary state or compatibility backends.
3. Preserve floating-point order, source splitting, limiters, regularization
   and all physical assumptions unless a separately authorized change
   justifies their alteration.
4. Physical-source evaluations remain cell-local and independent of grid
   indices. Mutable per-cell trial state belongs to the integrator/executor;
   accepted `q/qp` is published only after the complete Lie step succeeds.
5. A rejected step is fatal. There is no recoverable whole-step substrate
   rollback. Deferred geometry publication is a separate local transaction.
6. Retain the intentionally unsuffixed single-layer NAMELIST forms; they
   are not obsolete interfaces.

## Performance and regression gates

- Audit redundant `q -> qp` conversions and physical-property evaluations
  in hot loops before optimizing. Prefer preallocated workspaces and
  read-only per-layer configuration, without adding duplicate caches.
- Compare optimized and strict/debug builds, `TEST_STATE*`,
  `TEST_STATE_CONVERSION`, physical-source regressions, OpenMP
  determinism and representative complete examples.
- Measure performance with comparable builds/inputs before claiming a
  speedup. A compilation or lexical audit alone cannot certify numerical
  invariants or computational improvements.
- For the pending two-layer functional/scientific gates, use the
  [authoritative two-layer roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).

This is a current **contract**, not a list of unfinished steps from the
former state migration. For repository history and dated characterization
results, consult the immutable pre-cleanup snapshot above.
