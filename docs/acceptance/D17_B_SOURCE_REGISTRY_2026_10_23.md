# D17-B source metadata and legacy adapter — 23 October 2026

Status: **implemented, local GNU validation passed for phase B**. D17 as a whole
remains in progress. Baseline: `2faaa714` (committed D17-A inventory/contracts).
This implementation adds metadata beside the existing numerical path; it does
not yet route stage, boundary or mass-exchange evaluation through the registry.

Evidence: [machine-readable results](results_d17_b_source_registry_2026_10_23.json).
Reproduction: [standalone test README](../../TESTS/TEST_SOURCE_REGISTRY/README.md).
The [phase-A contract](D17_SOURCE_CONTRACT_2026_10_23.md) and its manifest retain
their original baseline evidence; their proposed IDs are now defined by the
catalog module. This report records the subsequent implementation separately.

## Implementation and compatibility

- `source_catalog_2d` defines all 34 stable IDs and their roles,
  representations, target blocks, real/complex capabilities, rheology membership
  and optional atomic-batch membership in one immutable table. Thirty entries
  are physical; four dependent assemblies are catalogued but rejected as
  independent registrations. Rheology family identity refers to the existing
  resolved model constants.
- `source_registry_2d` provides group/solver descriptors, separate stage/step
  registries, explicit identity/phase/order, and transactional initialization.
  Private storage is returned only through detached copies. Failure releases
  partial state; finalization is idempotent. A cached initialization flag avoids
  repeatedly traversing descriptor arrays.
- The legacy adapter validates the original velocity/energy/pressure group IDs
  against the existing compact equation partition. Shared Newton groups retain
  one jointly solved equation block. Shared analytic groups expand in the
  existing velocity → temperature → pressure order; split IDs start above four,
  with original IDs retained separately. Default EX stage ID four avoids the
  analytic group-one collision, while step ID one has a distinct namespace.
- Registration checks reject duplicate physical terms across cadences, unknown
  or dependent IDs, invalid methods/cadences/targets, incomplete or orphan solves,
  overlapping solve equations, unsorted/ambiguous execution order, mixed rheology
  families, partial analytic velocity bundles and incomplete/split catalogued
  atomic batches. The registry validates batch membership generically; it does
  not enumerate particular rheologies or physics terms. Active mass-exchange
  proposals remain in one step group.
- `source_configuration_2d` translates existing activation flags once after
  equation partition initialization. It records potentially active contributions;
  the existing physical kernels still apply wetness, packing, windows, reserve
  caps and per-cell source masks. Disabled rheology retains the zero-rate
  momentum block. An empty implicit partition has no solver descriptors.
- `simulation_context_type` owns and releases the registry. Its initialization
  reports invalid metadata before allocating the evolving state. The registry
  is rebuilt from existing input; no new namelist or restart payload is added.

Source rates, solver residual/Jacobian, tableaux, friction corrections, spatial
boundary handling and mass-exchange application are unchanged. The build graph
contains the new modules in dependency order, and `make distdir` distributes
both modules and the standalone test. The test is not added to GNU/ifx CI gates
in this phase; no GitHub Actions were launched.

## Local validation

GNU Fortran (Homebrew GCC 15.2.0) 15.2.0, macOS ARM64:

| Check | Evidence |
|---|---|
| D17 metadata, strict and optimized | 864 legacy configurations per profile; 14,823 assertions and 27 rejected invalid configurations per profile |
| Catalog-refactor follow-up, strict and optimized | 864 legacy configurations per profile; 14,857 assertions and 28 rejected invalid configurations per profile; logs in `build-local/d17-registry/refactor` |
| Standalone rheology suite, strict and optimized | 255/255 cell-local fixtures per profile; model names, source modules, analytic velocity, exact cooling/pressure and grouped DIRK checks pass |
| Configured full solver builds | optimized and strict-debug pass |
| Real parser plus short valid simulations | 32/32 checks per build; simulation-owned registry initialization exercised |
| Actual `n_RK=2` stage integration | 6/6 per build: quadratic, Coulomb stopping and lava, each at 1/4 threads |
| Source-mass suite | PASS: optimized, one thread; 10 mass-budget cases plus input rejection checks |
| Distribution/static checks | Source/test inclusion, manifest-ID match, shell syntax, dependencies and whitespace |

The catalog-refactor follow-up reran the focused strict and optimized registry
suite after replacing the per-ID descriptor switch with the catalog table and
making atomic-batch validation and test generation data-driven. The full-solver,
parser and integration results above are the original phase-B evidence; they
were not rerun for this metadata-only follow-up. The machine-readable JSON
retains the original phase-B run and its source hashes.

The strict standalone D17 profile uses `-O0 -g -Wall -Wextra -fcheck=all
-finit-real=snan`; optimized uses `-Ofast -funroll-all-loops`. Both enable
invalid/zero/overflow traps. Full solver build flags and raw-log hashes are
recorded in the JSON evidence. All generated objects/modules/fixtures remain
under `build-local/d17-registry` or test-owned temporary directories, outside
`src`. No tracked operational test/example input or numerical tolerance changed.

## Gate boundary and next steps

Accept phase B for this local scope: descriptors, initialization-time validation,
immutable registries, legacy normalization and simulation ownership are present.
This does **not** certify a new numerical integration path: the existing solver
still consumes `equation_partition_type` and its current ordered groups.

D17-C must add the analytic assembly bridge and remove physical law calls from
the numerical solver without changing state-refresh order or effective IM
increments. D17-D/E must connect term/solver filtering and explicit/boundary/step
paths while preserving coupled blocks and atomic mass proposals. D17-F must add
new native input syntax/presence tracking, remaining capability checks and full
numerical/portability acceptance. Unsupported partial bundles remain rejected;
this phase does not claim arbitrary term/method combinations are available.

Copy accessors and full registry validation are for initialization/diagnostics.
Future hot-path integration must bind/precompute the required metadata outside
cell/Newton loops; copying or validating the registry per cell would reopen the
D16 performance risk. The historical D16 ifx CI acceptance is not evidence for
these new D17 modules. No new Intel ifx/Linux/Docker acceptance is claimed.
