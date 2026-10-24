# D17 — Source-system refactoring: implementation roadmap and handoff log

**Repository:** [demichie/IMEX_SfloW2D_dev](https://github.com/demichie/IMEX_SfloW2D_dev)  
**Working branch:** `feature/d17-source-refactor` (branched from `master`)  
**Last updated:** 2026-10-24  
**Status:** IN PROGRESS — **not accepted**; GNU build **passed** after import fix, numerical regression tests not yet executed.

> **Read this file first when resuming work in another chat.** This document records the agreed *target architecture*, observed repository state, actual branch commits, CI evidence, and the next action. It is not a declaration that a change has passed compilation or numerical regression tests. Update this file in the same branch whenever D17 code changes or CI status changes.

## 1. Scope and authority

D17 reorganizes the configurable physical source terms and their time-integration groups, separating *source configuration*, *catalog metadata*, *source assembly*, *physical source laws*, *numerical solvers*, and *state maintenance*. Its primary acceptance criterion is preserving numerical and physical behavior with the legacy/default configuration.

The consolidated general roadmap already exists at [`docs/IMEX_SfloW2D_Development_Roadmap_WITH_LEGACY_LABELS_2026-10-15 (1).md`](./IMEX_SfloW2D_Development_Roadmap_WITH_LEGACY_LABELS_2026-10-15%20(1).md), particularly section D17; supplementary evidence is under [`docs/acceptance/`](./acceptance/). **Warning:** parts of that older roadmap describe pre-revision architecture (per-contribution IDs and multiple ordered Newton solvers). When they conflict with the explicit decisions below, the decisions below govern this branch. Do not silently reinterpret the older roadmap as the current target.

**Not yet present in this branch as a tracked source of truth:** the latest fully elaborated standalone D17 specification discussed in the planning conversation (D17.0–D17.21). Before the large architectural replacement, add or reconcile that complete specification without dropping detailed implementation, test and acceptance steps. This document is a self-contained implementation handoff and progress register, **not a substitute for the full technical specification**.

## 2. Architectural decisions — binding constraints

1. **One catalog entry per independently configurable physical source/model**, not one per contribution to its constitutive law. For example, `rheology_voellmy` owns both Coulomb and turbulent drag. Deposition incorporates the associated carrier deposition loss; do not keep its loss as a separate schedulable source.
2. **Minimal source catalog:** canonical `name`, `target` equation block (momentum, temperature, pore pressure as applicable), and `has_explicit`, `has_newton`, `has_analytic`. **Do not retain permanent numeric source IDs** or catalog metadata such as role, family, model_id, bundle_id, real_rate, complex_rate, has_semi_implicit. Numeric array positions, equation indexes, group/order identifiers, and method constants are permissible where numerically meaningful. Source identity is a name throughout.
3. Methods: `EXPLICIT`, `NEWTON`, `ANALYTIC`. Each physical source law exposes **four** entry points: explicit, semi-implicit, implicit, analytic. `NEWTON` uses both semi-implicit and implicit callbacks; semi-implicit may intentionally return zero for supported models. Unsupported method combinations must fail **at input initialization**; unsupported direct calls must fail clearly, not silently return zero. Keep the method mechanism extensible for future integration strategies, without implementing ML now.
4. **Static Fortran dispatch** in `rheology_laws/rheology_2d.f90` by canonical name with straightforward `SELECT CASE`; no factories, procedure-pointer registry, generated dispatcher, or extra generic layers. For a new conventional rheology, changes should be restricted to its law module, catalog entry, dispatcher branch, input parameters, `src/Makefile.am`, and tests; avoid changing generic source configuration, registry, assembly, Newton solver, integrator or state handling.
5. **Exactly one coupled Newton solve per NEWTON group.** Preserve the generic solver's equation-index lists and convergence/line-search machinery; eliminate the current repeated solver-descriptor sweep. **One analytic solve per ANALYTIC group**; use distinct groups for separate analytic solves.
6. **Momentum direction correction**: move `correct_friction_momentum` out of rheology into `src/state/state_2d.f90`, rename `correct_momentum_direction`. Apply after the semi-implicit state update **before** forming the Newton residual's semi-implicit term, and once after Newton convergence relative to the pre-Newton momentum state, only when momentum is included.
7. **State validity and dry-cell policies** belong in `src/state/state_2d.f90`, with correct call order from `time_integration_2d.f90`; avoid RK-specific naming for generic state actions. Keep conversion formulae in `state_conversion_2d.f90`. Include dry-cell primitive initialization in reconstruction halos.
8. Rename `src/physics/equation_terms_2d.f90` to `src/physics/source_terms_2d.f90` and provide the four general assemblers `eval_expl_terms`, `eval_semi_implicit_terms`, `eval_implicit_terms`, `eval_analytic_terms`. Keep constitutive formulae in the existing/specific laws modules and leave dispatchers free of physics, limits and numerical algorithms.
9. **Lava thermal cooling is independently configurable** from lava rheology. Add an activation flag in `TEMPERATURE_PARAMETERS`. If unspecified, preserve legacy default lava=>cooling on; explicit on/off overrides it for lava and nonlava. Retain existing cooling formulae and genuine rheology dependencies elsewhere.
10. **Mass exchange** retains its atomic inventory/bed/state update, shared limiter, and legacy late `STEP/EXPLICIT` cadence while its assembly is integrated into the common source interface. Proposed home: `src/physics/mass_exchange/{deposition_2d,erosion_2d,liquid_loss_2d,darcy_loss_2d,entrainment_2d,vaporization_2d,mass_exchange_laws_2d}.f90`. Do not enable unverified stage scheduling or Newton/analytic mass exchange.
11. Preserve the existing special treatment of boundary injection and its effective routing; no architecture work on that pathway is included in this change.
12. No additional registries, role/family descriptors, semi-implicit catalog capabilities, or other abstractions outside this plan. Keep defaults and accepted science intact.

## 3. Observed repository baseline at D17 branch creation

Branch created from commit `69b62511224e03dca30fd6a83b73ecafe3568bfa`. At that point earlier D17-related groundwork existed on `master`: source catalog with **34 old entries**, registry and native/legacy normalization, group-aware rate routing, a separately scheduled post-RK mass exchange module, and D17 acceptance/test artifacts. The existing architecture **does not yet satisfy** the revised minimal, name-only model. Older tests still assert legacy numeric IDs and original source count. These pre-branch commits must not be confused with work completed *on* the D17 branch.

Main files to coordinate:
- `src/numerics/source_catalog_2d.f90`
- `src/numerics/source_registry_2d.f90`
- `src/numerics/source_input_2d.f90`, `src/physics/source_configuration_2d.f90`
- `src/physics/equation_terms_2d.f90` (rename planned)
- `src/physics/rheology_laws/rheology_2d.f90`, per-model rheology modules
- `src/numerics/nonlinear_solver_2d.f90`, `src/numerics/time_integration_2d.f90`
- `src/numerics/mass_exchange_step_2d.f90`, `src/physics/mass_exchange_laws_2d.f90`
- `src/state/state_2d.f90`, `src/Makefile.am`, input parser, relevant `TESTS/`.

## 4. Actual commits on `feature/d17-source-refactor`

The following commits are **published changes**, not necessarily compiled or regression-validated:

| SHA | Change | Verification |
|---|---|---|
| [`fcb3d373`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/fcb3d373056d95b84bd0adcb977b09fe4599c9f2) | Added catalog lookup by source name, name uniqueness check, and `group_contains_name`. | Not independently tested |
| [`bd2d1cbc`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/bd2d1cbc58fd2c8040e5785ee60f300038ccdd85) | Materialized `term_names` alongside old `term_ids` in groups and solvers. | Transitional, must be simplified |
| [`922cb5e2`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/922cb5e2dc363fb5e086d152f6a5d293bf19c875) | Preserved direct validation compatibility for legacy group descriptors. | Not independently tested |
| [`d4ffc17a`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/d4ffc17a7dc2e3977cd3ab225b92edd0ed825c6c) | Added name lookup/group membership tests. | Tests added, not executed |
| [`27758c8f`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/27758c8fd8f159c884de2b91c3285f2f8fcc8947) | Added three method capabilities to the old catalog and used them in group checks. | Capabilities reproduce previous eligibility; not verified against every physical interface |
| [`bd3cd731`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/bd3cd731925d7d5b2ba34872f596d0a81a4e43cb) | Added capability-matrix tests. | Not executed |
| [`533cbe9d`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/533cbe9d8e5f769affbdb3a5f6bdbaafb52760f3) | Name-based descriptor resolution in registered group validation. | Not independently tested |
| [`447ce9a5`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/447ce9a57cca335c187c3cdf2cbd9bf7c40c441a) | Name/descriptor capability round-trip tests. | Not executed |
| [`66ea646e`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/66ea646e6e94dccf89e7330170348ac5991677aa) | Name-based selection for four cell-local explicit sources. | Not independently tested |
| [`9c908505`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/9c908505805391e09a2fa4baca8050d1f025b1d5) | Name-based Newton group membership and role metadata lookup in source assembly. | Initially failed compile due to missing imported symbol |
| [`fa1820a5`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/fa1820a53467c5db6fa846ca1b29bcdb1d517789) | Added `.github/workflows/d17-compile.yml`: GNU build-only CI, full log artifact even on failure. | First run failed as below |
| [`e3a9f743`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/e3a9f743faa4683c1804c4b81fdad08131d71078) | Added missing `source_term_descriptor_by_name` import in `equation_terms_2d.f90`. | **GNU compilation succeeded: run 40084588414** |

**Important:** The multiple temporary `term_ids`/`term_names` dual paths in these commits are a **migration scaffold only**. Do not extend this scaffold indefinitely: the agreed destination is a **single name-only representation**. Existing old source sub-contributions and over-rich catalog fields must be removed in coordinated commits.

## 5. Compilation evidence and actionable diagnostics

Workflow: [`D17 compile diagnostics`](../.github/workflows/d17-compile.yml). Its trigger is every push to `feature/d17-source-refactor`, plus manual dispatch; it runs `autoreconf -fi`, `./configure --with-netcdf=/usr`, and `make -j2`, **without** tests. On any outcome it attempts to upload `environment.log`, `autoreconf.log`, `configure.log`, `compile.log` and `config.log`.

**First recorded run:** [run #1, ID 40084188010](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084188010), against `fa1820a5`, **FAILURE**.
- Ubuntu 24.04 / GNU Fortran 13.3.0, netCDF-Fortran 4.5.4.
- Compiler and dependencies installed; Autotools generation and `configure` successful.
- Compilation halted in `src/physics/equation_terms_2d.f90:737`, at `source_term_descriptor_by_name(solver%term_names(i))`: missing import in `USE source_catalog_2d, ONLY`. Diagnostic included both `has no IMPLICIT type` and a PURE-procedure purity error (latter is a consequence of missing interface). No evidence of independent physics/solver errors in this run.
- Fixed by `e3a9f743` with the missing import. **Follow-up CI run on e3a9f743: SUCCESS** — [run 40084588414](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084588414), compile-only, no test execution.

**Next CI action:** confirm build status after the documentation commit; then run the existing source-registry test suite separately. Compile-only SUCCESS does not certify numerical equivalence. If a later build fails, analyze the first error in its `compile.log`.

## 6. Work remaining — roadmap and sequencing

Status legend: **NOT STARTED**, **PARTIAL**, **PENDING CI**, **VERIFIED** (only when documented evidence exists).

| Work package | Status | Required deliverable / gate |
|---|---|---|
| D17 specification consolidation | **PARTIAL** | Commit/reconcile the latest full detailed D17.0–D17.21 design and implementation plan with this handoff. |
| GNU build gate | **VERIFIED for e3a9f743** | [Successful compile-only run 40084588414](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084588414); further commits still require their own check. |
| Catalog normalization | **PARTIAL** | One name-only descriptor per independently configurable source; only name, target, three capabilities. Remove permanent source IDs, old roles, families, bundles, etc. |
| Reology source consolidation | **NOT STARTED** | One source per rheology, internal force components evaluated by model; preserve model-specific formulae. |
| Name-only groups and solvers | **PARTIAL** | Remove dual `term_ids`/`term_names`; eliminate numeric source ID conversion, update input/registry/solvers/tests together. |
| Method contracts and static dispatch | **PARTIAL** | Exact four law entrypoints, explicit/Newton/analytic capability validation, deterministic static `SELECT CASE` by name; no dynamic registration. |
| Newton group semantics | **NOT STARTED** | One coupled Newton system per group with supported equation list; preserve complex-step/line-search path. |
| Momentum correction relocation | **NOT STARTED** | Move and rename to `state_2d`, enforce pre- and post-Newton momentum direction corrections. |
| Generic state/dry-cell consolidation | **NOT STARTED** | Move state checks, bounds, roundoff corrections and halo dry-cell initialization to state module. |
| Source assembler/module migration | **NOT STARTED** | `source_terms_2d` filename and four assemblies; physical formulas only in respective laws. |
| Independent lava thermal cooling | **NOT STARTED** | New `TEMPERATURE_PARAMETERS` flag; absent flag maintains lava auto-on; four combinations tested. |
| Mass exchange modularization | **NOT STARTED** | Separate physical mass laws; one deposition plus carrier loss; shared limiter; post-RK atomic STEP/EXPLICIT. |
| New source demonstration | **NOT STARTED** | Show that adding standard rheology touches only model module, catalog, static dispatch, parameters/build/tests. |
| Source-registry strict/optimized regression | **FAILED / ROOT CAUSE PENDING LOG** | [Run 40084977002](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084977002) of the existing registry script; inspect result and update this row. |\n| End-to-end regression matrix | **NOT STARTED** | Build, source registry, input defaults, coupled Newton, real/complex/complex-step Jacobians, RK, dry limits, mass budget, cooling permutations, legacy run comparisons. |
| Final D17 acceptance | **NOT STARTED** | Complete detailed D17.21 evidence matrix and explicit review; do not declare done early. |

### Recommended next steps

1. **CI compile is green after the import fix.** Next run existing source-registry tests and compare legacy behavior before larger changes.
2. **Reconcile/commit the latest complete D17 specification** so future workers have one detailed authoritative implementation contract.
3. Prepare one **coordinated name-only change set** covering catalog, registry, source configuration, source assemblers, mass-exchange selectors and tests. Replace old per-contribution rheology items with one entry per configured source; avoid another layer of permanent adapters.
4. Implement static rheology dispatch and the four explicit/newton/analytic-compatible interfaces with one coupled Newton solve per group.
5. Implement state/momentum corrections, reorganized assembly and mass laws, independent cooling and regression gates in the D17 order.

## 7. Focused source-registry regression CI

Dedicated workflow: [`d17-source-registry.yml`](../.github/workflows/d17-source-registry.yml), added by [commit `788bbc70`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/788bbc704eec6c021ba16573710ee7e8f8245ba0). It runs the **existing** `TESTS/TEST_SOURCE_REGISTRY/run_test.sh` with GNU Fortran in both `strict` and `optimized` profiles; no new test logic or scientific semantics were introduced. Logs and individual profile outputs are uploaded as an artifact even on failure. The workflow triggers on D17 branch pushes or manually.

First run: [40084977002](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084977002), triggered by `788bbc70`. **Outcome: FAILURE**; checkout and compiler setup succeeded, the `Run existing D17 registry regression` step failed, and artifact upload succeeded. Root cause NOT established; inspect the artifact before changing source code. Check the live result, the run's job steps and `d17-registry-logs-*` artifact. If failing, investigate test output and record the first substantive cause; compilation-only CI success is not equivalent to passing source registry tests.

## 8. Instructions to the next chat/maintainer

- Read this handoff, the general roadmap's D17 section, and **the latest reconciled detailed specification once added**; inspect the current **branch HEAD**, never assume the last SHA here is still current.
- Review actual commits and GitHub Actions runs to know what changed and what passed. “Committed,” “compiles,” “unit-tested” and “numerically equivalent” are distinct statuses.
- Work **only** on `feature/d17-source-refactor`; keep `master` unchanged. Do not create a PR or merge without authorization.
- Preserve the user's agreed architecture, default physical behavior, dry-state semantics, atomic mass exchange, model-specific laws and coupling.
- After each meaningful commit, update **this file** with SHA, affected modules, CI URL/result, tests run and outstanding work, in the same branch; never mark an unrun test as passed.
- Avoid superficial one-function commits that accumulate temporary adapters without converging on the final design. Favor coherent, reviewable changes supported by CI.
