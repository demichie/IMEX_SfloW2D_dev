# Adding a rheology model to IMEX-SfloW2D

This guide describes the **current** source-group architecture on the development branch.
It is a developer tutorial, not an already-implemented constitutive law. The illustrative
`linear_drag` model below is **not** registered in the repository. Do not copy only its
physics module and expect the solver to accept the input.

The guiding separation is:

- `src/physics/rheology_laws/` owns constitutive formulas and their REAL/COMPLEX
  implementations.
- `src/numerics/source_catalog_2d.f90` describes source identities, available
  integration methods, phase permissions, and Newton equation footprints.
- The existing source configuration, registry, Lie executor, and nonlinear solver
  normally **do not** need model-specific cases.
- `src/io/inpout_2d.f90` owns rheology NAMELIST parsing and coefficient validation.
- `src/Makefile.am` owns source lists and explicit Fortran module dependencies.

Do not add global coefficient mirrors, per-cell configuration copies, alternate source
registries, compatibility dispatchers, or model-specific timestepping code.

## 1. Start with the physical law

As a deliberately simple example, introduce a linear momentum relaxation model.
For a mixture with thickness `h`, density `rho_m`, and horizontal velocity
`u_h = (u,v)`, define the conservative horizontal momentum

\[
\mathbf{m} = \rho_m h\,\mathbf{u}_h
\]

and the **cell-local rheological source**

\[
\left.\frac{d\mathbf{m}}{dt}\right|_{\mathrm{rheo}}
  = -\gamma\,\mathbf{m}, \qquad \gamma \ge 0,
\]

where `gamma` has units of `s^-1`. This is not a momentum flux and does
not change mass, species inventories, or thermal energy. It is an example of
linear drag, **not** a claim that such a law is appropriate for a particular flow.

For fixed `gamma`, density, and thickness, the exact source-only solution is

\[
\mathbf{u}_h(t+\Delta t)
  = e^{-\gamma\Delta t}\,\mathbf{u}_h(t).
\]

An explicit Euler substep instead gives
`m_new = (1 - gamma*dt)*m_old`, and can reverse the direction if
`gamma*dt > 1`. Therefore an explicit implementation must have a defined
stopping/stability policy; it is not automatically equivalent to the analytic
solution.

**State and numerical invariants:**

1. For `gamma >= 0`, the source dissipates horizontal kinetic energy;
   it must not inject horizontal momentum at rest.
2. Mass and composition stay unchanged by this source.
3. Exactly dry cells and `dt = 0` must be safe.
4. All active coefficients must be finite and dimensionally consistent.
5. Any REAL/COMPLEX pair must implement the same mathematical source.
6. Layer-specific parameters come from the supplied `rheology_config_type`.

Decide which methods are **actually implemented and tested** before advertising
them: `EXPLICIT`, `NEWTON`, and/or `ANALYTIC`. The smallest introductory
change supports `EXPLICIT` only. Other treatments are optional extensions,
not stubs to advertise prematurely.

## 2. Files to create and modify

| File | Action | Purpose |
| --- | --- | --- |
| `src/physics/rheology_laws/rheology_linear_drag_2d.f90` | **Create** | Implement the law and the supported treatment procedures. |
| `src/physics/rheology_laws/rheology_parameters_2d.f90` | Modify | Add the one canonical per-layer coefficient. |
| `src/physics/rheology_laws/rheology_2d.f90` | Modify | Import and dispatch the law in supported methods. |
| `src/numerics/source_catalog_2d.f90` | Modify | Register a source, its ID and its method/phase capabilities. |
| `src/io/inpout_2d.f90` | Modify | Validate the new coefficient and model; retain the existing NAMELIST mechanism. |
| `src/Makefile.am` | Modify | Add the compilation unit and its dependencies. |
| `TESTS/TEST_RHEOLOGY_CELL_LOCAL/` | Extend tests | Check physical rates, naming, parser validation, and integration. |
| `TESTS/TEST_SOURCE_REGISTRY/` | Extend tests | Check source groups and method/phase validation. |
| `EXAMPLES/` or a suitable existing test fixture | Optional | Provide an actual input and a representative simulation. |

Do **not** create a new source registry or new `numerics` executor. No changes
to `source_configuration_2d.f90`, `source_registry_2d.f90`, or
`nonlinear_solver_2d.f90` should be needed for an ordinary new law supported
by existing methods. Inspect exceptions rather than assuming none exist:
`lie_source_executor_2d.f90` currently contains a named exception for isolated
explicit `rheology_mui` momentum reversal. A new explicit law does not
automatically inherit that exception.

## 3. Add the coefficient to the canonical configuration

Edit `src/physics/rheology_laws/rheology_parameters_2d.f90`, adding a member
to `rheology_config_type`, for example:

```fortran
REAL(wp) :: linear_drag_gamma=-1.0_wp  !< Linear drag rate [s^-1].
```

The negative default is an **unset sentinel**, not a valid physical coefficient.
Both `rheology_l1` and `rheology_l2` already use this type: do not introduce
a second variable or new per-layer storage. Existing NAMELIST groups
(`rheology_parameters` for single layer, and the corresponding
`rheology_l1_parameters`/`rheology_l2_parameters`) bind the canonical
objects. Check their read/write paths when adding a member.

## 4. Create the physical module (EXPLICIT-only starting point)

Create `src/physics/rheology_laws/rheology_linear_drag_2d.f90`:

```fortran
!> @brief Linear conservative-momentum relaxation for one physical cell.
MODULE rheology_linear_drag_2d
    USE parameters_2d, ONLY : wp
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE state_types, ONLY : physical_state_real
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_linear_drag_explicit

CONTAINS

    !> @brief Evaluate the explicit conservative momentum source.
    !> @param[in] config Canonical configuration of the selected layer.
    !> @param[in] state Cell-local physical state.
    !> @param[out] rate Horizontal conservative momentum rates.
    SUBROUTINE eval_rheology_linear_drag_explicit(config,state,rate)
        TYPE(rheology_config_type), INTENT(IN) :: config
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(OUT) :: rate(2)

        rate=-config%linear_drag_gamma*state%rho_m*state%h*[state%u,state%v]
    END SUBROUTINE eval_rheology_linear_drag_explicit

END MODULE rheology_linear_drag_2d
```

This module returns a *rate*, not an updated velocity. The source assembler
in `src/physics/source_terms_2d.f90` and the Lie executor own rate assembly,
admissibility, and state commit. The module has no grid indices, shared mutable
scratch arrays, I/O, or `NAMELIST` parsing. The example expects coefficients
to have been validated before entering the physical kernel.

The derivation and the exact source-only update above provide two independent
checks: `rate == -gamma*m`, and, for fixed state, `u_after == exp(-gamma*dt)*u_before`.
They are *tests to implement*, not evidence that the illustrative module is
currently covered by CI.

## 5. Register the source in the catalog

Edit `src/numerics/source_catalog_2d.f90`. Register the canonical source name
`rheology_linear_drag` and the user-facing model name `linear_drag`.
For the EXPLICIT-only implementation, the descriptor should have:

```fortran
source_term_type('rheology_linear_drag',.TRUE.,.FALSE.,.FALSE., &
    default_method=SOURCE_METHOD_EXPLICIT, &
    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.)
```

The phase permissions above are **illustrative**. Select them from the model's
actual physical domain of validity; do not copy permissions without justification.

**Important catalog invariant:** the current catalog starts with a contiguous
block of rheology entries, numbered `RHEOLOGY_ID_VOELLMY` through
`RHEOLOGY_ID_MUI_VOELLMY`. `rheology_source_id` searches that initial block.
Adding a new rheology **after** non-rheological sources will make this resolver
fail unless you redesign the lookup. The minimally invasive approach today is:

1. Insert the new descriptor immediately after the existing rheology entries.
2. Assign `RHEOLOGY_ID_LINEAR_DRAG` to the corresponding catalog position.
3. Update `SOURCE_TERM_COUNT` and **all shifted fixed
   `SOURCE_ID_*` constants**, including mass-exchange and interlayer entries.
4. Update the bound in `rheology_source_id` so the new entry is included.
5. Audit every index range, source bit mask, fixed-ID assumption, and test
   involving the catalog order. IDs are internal, not stable external ABI.
6. Keep the catalog within the available integer bit-mask capacity; increasing
   source count is not safe merely because `SOURCE_TERM_COUNT` was increased.

The generic function `canonical_rheology_source_name` already derives
`linear_drag` from a canonical name prefixed by `rheology_`. The registry
then caches the selected rheology ID. Neither needs a new model-specific switch.

## 6. Connect the physical dispatcher

In `src/physics/rheology_laws/rheology_2d.f90`:

- Import `RHEOLOGY_ID_LINEAR_DRAG` from the catalog.
- Import `eval_rheology_linear_drag_explicit` from the new module.
- Add a `CASE(RHEOLOGY_ID_LINEAR_DRAG)` in `eval_rheology_explicit`:

```fortran
CASE(RHEOLOGY_ID_LINEAR_DRAG)
    CALL eval_rheology_linear_drag_explicit(config,state,process_rate)
```

For EXPLICIT-only registration, **do not advertise** unsupported `NEWTON`
or `ANALYTIC` capabilities in the source catalog. Keep the shared
`EXPLICIT`, `SEMI_IMPLICIT`, `IMPLICIT` and `ANALYTIC` law interfaces:
an intentionally absent contribution may return zero; an unsupported
treatment must report its established explicit error. Inspect the
exhaustive `SELECT CASE` statements in `eval_rheology_semi_implicit_id`,
`eval_rheology_implicit_id` and `eval_rheology_implicit_real_id` and
add the new ID to the appropriate zero/unsupported cases. Do not remove
these common interfaces as an apparent code-cleanup opportunity.

If implementing `NEWTON` as well, provide the physically correct real and
complex implicit evaluations (or a justified split between semi-implicit and
implicit parts); add the model to `source_newton_affected_equations` and
`rheology_name_has_role` in the catalog, and ensure the dispatcher calls the
correct REAL and COMPLEX procedures. The affected equations for simple momentum
drag are `idx_uEqn` and `idx_vEqn`. Check that Newton residual/Jacobian
and complex-step tests agree. Do not introduce an independent approximation
solely for the complex path.

If implementing `ANALYTIC` as well, provide an exact source-flow procedure and
a matching case in `dispatch_analytic_rheology`. For the simple linear law,
use the exponential update above. Mark `has_analytic=.TRUE.` **only after**
its implementation and tests exist. A model that changes mass or composition
during the source step needs a separate analysis; this example does not.

## 7. Extend input validation

Edit `src/io/inpout_2d.f90`. The current file owns the NAMELIST declarations,
single/two-layer input handling, and the model-specific coefficient checks.

For the example:

1. Include `linear_drag_gamma` in the existing rheology coefficient access
   or explicit canonical configuration validation.
2. Add `linear_drag` to the active-model finite-coefficient check
   (using the repository's `finite_input_value`).
3. Reject negative values and the unset sentinel; zero is valid and means
   no drag.
4. Update model-specific validation/diagnostics and any supported-model
   strings that are maintained by the parser.
5. Confirm `EXPLICIT` is accepted, while `NEWTON` and `ANALYTIC`
   are rejected through the source catalog when not implemented.
6. Test both the short single-layer NAMELIST and L1/L2 configuration rules.
   The short unsuffixed NAMELIST remains a deliberate single-layer interface,
   not a transitional alias to remove.

A minimal **illustrative** single-layer input fragment is:

```fortran
&RHEOLOGY_PARAMETERS
  rheology%rheology_model='linear_drag',
  rheology%rheology_method='EXPLICIT',
  rheology%rheology_group=1,
  rheology%linear_drag_gamma=0.5,
/
```

This is **not** a complete runnable input file. Follow an existing working
example for phase flags, material parameters, domain, initial condition, output,
and other mandatory settings; the source should be assigned a valid Lie group.
Do not change other physical assumptions merely to exercise the new law.

## 8. Update the build dependency graph

Edit `src/Makefile.am` to add the new `.f90` to the solver source list and to
declare all `USE`-module prerequisites for its object. Also make
`rheology_2d.$(OBJEXT)` depend on the new law's object. Inspect the nearby
`rheology_quadratic_2d` rules as a structural example.

Run `autoreconf -fi` when regenerating the Autotools build. A successful
single-threaded compilation is not proof of a correct parallel dependency
graph: test `make -j4` from a clean build.

## 9. Add focused tests before integration tests

Extend the existing suites rather than creating a second test framework:

1. `TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_model_names.f90`: canonical model name,
   unique registration, method availability and phase permissions.
2. A cell-physics test under `TESTS/TEST_RHEOLOGY_CELL_LOCAL/`: zero velocity,
   `gamma=0`, two nonzero velocity components, exact sign, dimensions,
   dry state, and the analytic rate `-gamma*m`.
3. `TESTS/TEST_RHEOLOGY_CELL_LOCAL/check_inputs.py`: valid input and rejection
   of missing, negative, NaN and infinite `gamma`; reject unsupported methods.
4. `TESTS/TEST_SOURCE_REGISTRY/test_source_registry.f90`: correct EXPLICIT Lie
   group and correct rejection of unsupported methods and phase combinations.
5. If adding a Newton implementation: verify both REAL/COMPLEX evaluations,
   complex-step Jacobians and the coupled stage residual.
6. If adding an analytic implementation: test exact exponential decay,
   zero-step identity, non-increasing speed and the semigroup property
   `flow(dt1+dt2) == flow(dt2,flow(dt1))`.
7. Run a representative single-layer simulation, check conservation of mass,
   physical admissibility, wet/dry front behavior, and OpenMP reproducibility.
   For a two-layer-enabled configuration, verify independent coefficients and
   layer-specific results.

An explicit friction law can generate a reversing Euler candidate when
`gamma*dt>1`. The current Lie executor includes a **mu(I)-specific**
isolated-source reversal projection, not a generic automatic guarantee for a
new law. Test the actual explicit update and obtain an agreed numerical
treatment for overshoots **before** declaring the new model production-ready.
Do not silently alter the spatial operator, physical law, or timestep policy.

Suggested validation entry points (adjust paths to the actual build layout):

```sh
python3 tools/docs/check_fortran_style.py src/
python3 TESTS/TEST_SOURCE_DOCUMENTATION/check_documentation.py
sh TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh
sh TESTS/TEST_SOURCE_REGISTRY/run_test.sh
# Then clean optimized and strict-debug builds, full make check,
# representative examples, and relevant GitHub Actions.
```

Record numerical evidence, not just successful compilation. Do not claim CI
validation based on a run for an earlier commit.

## 10. Review checklist

- [ ] The physical equation, units, active phase range and dry-state limits are documented.
- [ ] There is one canonical coefficient per configured layer, with strict input validation.
- [ ] Physical source procedures use immutable explicit configuration arguments.
- [ ] All advertised treatments actually exist; REAL/COMPLEX methods agree where required.
- [ ] Catalog IDs, bit masks, count, lookup bounds and Newton footprints are consistent.
- [ ] Existing input, source grouping, Lie execution and restart semantics are preserved.
- [ ] `src/Makefile.am` includes the new source and every new module dependency.
- [ ] Focused tests, optimized/strict builds, numerical regressions and OpenMP checks pass.
- [ ] No redundant dispatcher, source registry, mutable global copy or compatibility layer was added.

## Existing implementations to study

- `src/physics/rheology_laws/rheology_quadratic_2d.f90`:
  existing physical kernels, REAL/COMPLEX sharing, and analytic treatment.
- `src/physics/rheology_laws/rheology_froude_2d.f90`:
  a law with EXPLICIT and ANALYTIC capabilities but no Newton capability.
- `src/physics/rheology_laws/rheology_mui_2d.f90`:
  a more complex stress-dependent law; do not use its regularizations by default.
- `docs/FORTRAN_STYLE.md`: mandatory lexical indentation, 120-column source
  width, Doxygen documentation and verification commands.

Every model has its own physical assumptions; this guide describes **integration
interfaces**, not a license to copy another law's numerical regularization.
