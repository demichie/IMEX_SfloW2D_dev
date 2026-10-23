# D17-B source registry and legacy normalization

Run from the repository root with GNU Fortran:

```sh
D17_WORK_DIR="$PWD/build-local/d17-registry/focused" \
  sh TESTS/TEST_SOURCE_REGISTRY/run_test.sh
```

The standalone runner retains objects, modules, executables and `results.log`
under the supplied directory (a retained `/tmp/imex-d17-registry.*` directory
by default). It runs strict (`-O0`, runtime checks, initialized signaling NaNs)
and optimized (`-Ofast`) profiles with invalid/zero/overflow traps. `FC` selects
a GNU-compatible compiler. No solver input, tableau or GitHub workflow changes.

Each profile checks all 864 combinations of the three legacy group assignments,
the two supported methods and optional energy/pore targets. Assertions preserve
joint Newton equation blocks, analytic velocity/temperature/pressure order,
public legacy IDs, collision-free split IDs, separate stage/step namespaces,
source activation, immutable copy access and idempotent lifecycle. It also tests
multiple independent Newton blocks and rejects invalid groups, solver references,
capabilities, dependent children, partial analytic bundles and incomplete
catalogued atomic batches. Rheology families and atomic batches are enumerated
from the catalog, so adding a family or batch does not require extending a
hand-maintained test list. The executable prints the assertion and
rejected-configuration counts.

`source_catalog_2d` defines the 34 D17 IDs/capabilities; dependent assembly IDs
are catalogued but cannot be scheduled. `source_registry_2d` owns validated group
and solver copies. `source_configuration_2d` translates existing configuration
once, and `simulation_context_type` owns/finalizes the registry. Explicit group
4 records cell rates and the boundary-state hook; step group 1 retains the shared
bounded mass-exchange transaction. Per-cell state/window/packing guards remain
in the existing physics. Explicit footprints conservatively retain the legacy
full equation vector, while implicit solve descriptors use exact target blocks.

The registry now drives analytic and explicit stage routing, ordered implicit
solves, and validation/execution of the atomic post-RK mass-exchange batch.
The optional `SOURCE_GROUPS` namelist is parsed and validated; the input-path
integration test covers a representative native schedule, not every valid or
invalid physical configuration. Boundary-source evaluation continues through
the equation bridge. The registry is therefore more than metadata, but the
remaining D17.20 matrix and cross-platform acceptance are still required.
The committed phase-A manifest remains historical evidence; phase-B and later
reports record the implementation and their respective validation scopes.

## Registering a rheology in D17-B

The phase-B registry is metadata-only; the existing model ID, law module,
build list and `rheology_2d` dispatch still need their normal updates. Register
each source contribution by adding its stable ID and one row to
`SOURCE_CATALOG`, including its `model_id` and capabilities. The existing
`rheology_source_terms` lookup derives the family from those rows. If multiple
contributions must share one kernel evaluation, assign them the same nonzero
`bundle_id`; the generic registry validator and this test discover that batch
automatically. No per-rheology cases are needed in `source_registry_2d`,
`source_configuration_2d` or the generated legacy-combination matrix. Add
physics-specific numerical tests for the new law as usual.
