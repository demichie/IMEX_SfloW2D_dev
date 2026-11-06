# Debug compiler warnings: declaration-only cleanup (2026-11-06)

## Scope and baseline

Baseline: user-supplied gfortran strict-debug build log using
`-O0 -g -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
The recorded compilation completed, with **198 warnings**:

| Category | Before | Action |
| --- | ---: | --- |
| Unused dummy arguments | 105 | Preserve current public/internal call signatures; assess as API maintenance separately. Most originate from treatment entry points that intentionally reject unsupported methods. |
| Unused `USE, ONLY` imports | 49 | Remove precisely the unused names from their host procedures/program. |
| Unused local declarations | 29 | Remove precisely the declarations; do not change executable expressions or physics. |
| `INTENT(OUT)` arguments not set | 15 | Preserve intentionally failing treatment interfaces; the bodies invoke `fatal_error`. |

### Files changed

| File | Unused imports removed | Local names removed |
| --- | ---: | ---: |
| `src/physics/mass_injection_terms/mass_bottom_circular_injection_2d.f90` | 12 | 0 |
| `src/physics/mass_injection_terms/mass_bottom_fissural_injection_2d.f90` | 15 | 0 |
| `src/physics/source_terms_2d.f90` | 19 | 1 (`pressure_rate`) |
| `src/numerics/lie_source_executor_2d.f90` | 2 | 0 |
| `src/IMEX_SfloW2D.f90` | 1 | 0 |
| `src/state/state_conversions.f90` | 0 | 4 (`rho_c`, `inv_rho_c`, `alphas`, `alphag` in the complex conversion) |
| `src/spatial/reconstruction_2d.f90` | 0 | 23 (unused `reconstruct_kt_eta_traces` temporaries) |
| `src/io/inpout_2d.f90` | 0 | 1 (`i_solid`) |
| **Total** | **49** | **29** |

The removed items were identified by the compiler's exact warning
locations. No procedure argument lists, physical formulas, time
integration logic, statements executable at runtime, `.inc` kernels
or model interfaces were changed.

## Warnings intentionally retained

Many unused arguments are required by the four physical source
evaluation interfaces: explicit, semi-implicit, implicit and analytic.
The unimplemented entry points call `fatal_error` and must remain
distinct from a successful zero-rate evaluation. For example,
`eval_pore_pressure_darcy_implicit` and
`eval_bottom_circular_injection_analytic` intentionally reject their
unsupported methods. `fatal_error` ends in `ERROR STOP 1`.
Forcing an output to zero solely to suppress `INTENT(OUT)` warnings
would obscure this contract and is not part of this change.

A smaller set of warnings concern unused arguments in implemented
paths (e.g. `apply_hp_reconstruction`, `lnsrch`,
`eval_pccu_interface`, `pore_gas_loss_complex`). Their signatures
were left intact pending call-site/API review, rather than modifying
a numerical interface opportunistically.

## Verification and remaining work

- Checked against the original source that, after excluding Fortran
  `USE` statements and local declarations, **all other source
  lines are identical** in the eight modified files. This verifies
  that this patch only changes declarations and imports, not
  executable numerical statements.
- Checked edited multiline `USE, ONLY` continuations for required
  commas and ampersands in column 120; any unrelated, pre-existing
  formatting issues were left unchanged.
- No floating-point operations or ordering were changed.
- **A fresh debug compile and numerical regressions have not been
  run** in this pass. Run the manual **Build / debug** Action to confirm
  the expected reduction from 198 warnings to approximately 120,
  assuming the compiler/settings match the supplied log. Treat any
  unexpected new diagnostics as blockers for further cleanup.
- After a successful compilation, relevant `TEST_STATE`,
  `TEST_SPATIAL_OPERATOR`, `TEST_MASS_EXCHANGE` and
  `TEST_SOURCE_MASS` regressions should be run manually.
- A later, separate API review may reduce unused-dummy warnings,
  but must preserve each source law's supported/unsupported
  treatment contract. Do **not** silence `-Wall` globally.
