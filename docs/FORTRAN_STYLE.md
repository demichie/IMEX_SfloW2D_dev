# Fortran source style (mandatory)

Applies to all production `src/**/*.f90` and to Fortran tests when edited.

## Indentation and layout

- Four ASCII spaces per lexical nesting level; **no TAB characters**, including in comments. Fortran keywords, logical constants, intrinsic logical operators and named relational operators use uppercase; identifier names retain their spelling and case.
- The indentation is **relative to the actual Fortran construct**: module-level USE/IMPLICIT/PUBLIC/PRIVATE/declarations are exactly 4 columns; module procedures and their corresponding END statements start at 4 columns; procedure declarations and top-level executable statements start at 8 columns. A multiple of four alone is not sufficient.
- For nested IF/ELSE IF/ELSE/END IF, DO/END DO, SELECT CASE/CASE/END SELECT, WHERE/ELSEWHERE/END WHERE, BLOCK and ASSOCIATE, each nested body adds exactly 4 spaces. CASE and CASE DEFAULT align with their SELECT, and statements inside each CASE are indented one additional level.
- Place a blank line before and after major IF, DO, SELECT, WHERE, BLOCK and ASSOCIATE constructs where adjacent executable statements would otherwise touch them. Keep CASE labels and ELSE branches attached to their associated construct; do not insert blank lines inside a continued statement or between a procedure's Doxygen header and its declaration.
- Apply indentation, blank-line grouping, width limits, argument packing and continuation alignment as **one coherent formatting pass**, and verify all these invariants together after every formatting change.
- Align `END IF`, `END DO`, `END SELECT`, `END BLOCK` and `END SUBROUTINE name` with their respective openers.
- Indent continued statements consistently with spaces; preserve the semantics of free-form `&` continuation and OpenMP/compiler directives.
- Separate declaration blocks, initialization, state conversion, source evaluation, admissibility checks, commit/rollback and cleanup with blank lines. Avoid stacked unrelated statements.
- Separate procedures visibly with a blank line and their own Doxygen header; use meaningful section comments rather than decorative banners inside procedures.
- No trailing whitespace. Pure formatting changes may normalize keyword/operator case but must preserve identifiers, literal values, the remaining executable tokens, floating-point expression order, runtime behavior and OpenMP directives.


## Line width and continuation policy

- Fortran executable code uses a maximum width of 120 columns.
- An existing trailing continuation marker is aligned at column 120.
- Pack as many complete arguments as will fit on each line; never split identifiers.
- Doxygen prose is wrapped at a maximum width of 100 columns.
- OpenMP directives and character literal continuations require syntax-aware handling.
- Pure layout changes must preserve all executable tokens apart from intentional keyword/operator case, and must preserve expression order.

## Mandatory Doxygen coverage

Document every module, derived type and significant field or variable (meaning, units, rank, ownership when relevant). Every FUNCTION and SUBROUTINE, including PRIVATE, PURE, ELEMENTAL, internal and type-bound procedures, must have an immediately preceding `!> @brief` header, `@details` when needed, and **one correctly directed** `@param[in]`, `@param[out]` or `@param[in,out]` for every dummy argument, including `this`; functions also require `@return`. Describe dimensions, units, preconditions, error behavior, side effects and thread-safety where applicable. Never invent generic documentation to satisfy a checker.

Comment substantial algorithmic blocks with the reason for ordering and key invariants. Do not narrate trivial assignments.

## Checks

The repository already provides `python3 TESTS/TEST_SOURCE_DOCUMENTATION/check_documentation.py` for Doxygen contract checks, including argument names/directions, derived-type `@brief` headers and Doxygen component comments;
it also supports `--check-against-git REF` to compare executable lines. The complementary `python3 tools/docs/check_fortran_style.py src/` checks four-space indentation, TABs and trailing whitespace. Run `python3 tools/docs/check_fortran_style.py --structural-report src/` to collect **non-failing, conservative** warnings for mismatched DO/IF/SELECT/WHERE block indentation. This heuristic intentionally does not rewrite source and may miss or misclassify continued statements; review its findings against the Fortran control-flow structure before editing. The original production sources received an initial whitespace-only normalization pass. The documentation suite invokes the style checker over `src/` on every run. Subsequent semantic reviews must still verify true lexical nesting, blank-line grouping and full Doxygen completeness; a multiples-of-four check alone does not prove these.

Formatting is not an opportunity to alter formulas or the Lie algorithm.

## Official Fortran formatter: Codee Formatter

[Codee Formatter](https://www.codee.com/codee-formatter/) is the official
formatter for Fortran source in this repository. Its configuration is maintained
in the root-level [`codee-format`](../codee-format) file. This guide remains
authoritative where formatter options cannot express the complete convention.
The earlier fprettify experiment has been retired.

### Automated formatting rules

- Four ASCII spaces for every lexical nesting level and continuation indent.
  The indentation must match the actual Fortran construct, not simply be a
  multiple of four.
- Maximum free-form width: 120 columns. Align trailing continuation `&`
  markers at column 120, subject to manual review for special constructs.
- **Uppercase Fortran keywords** (`MODULE`, `SUBROUTINE`, `IF`, `END IF`),
  logical constants (`.TRUE.`), logical operators (`.AND.`) and named
  relational operators (`.EQ.`). Preserve identifier names and their case.
- No more than **one consecutive blank line**; exactly **one blank line between
  procedures**; no extra blank lines at the beginning or end of a file.
  This is implemented with `ConsecutiveEmptyLines`.
- Remove trailing whitespace and retain a final newline. Preserve the current
  choices for declaration syntax, operator spacing and numerical notation.

### Using Codee Formatter

After installing Codee Formatter, inspect the resolved configuration from the
repository root:

```sh
codee format --config codee-format --dump-config
```

Confirm `IndentSize: 4`, `ColumnLimit: 120`,
`AlignAmpersandToColumnLimit: true`, `Casing.Keywords: Uppercase`,
`ConsecutiveEmptyLines.MaxToKeep: 1` and
`ConsecutiveEmptyLines.BetweenProcedures: 1`.

Preview the result on a copy instead of changing the repository immediately:

```sh
cp src/numerics/source_registry_2d.f90 /tmp/source_registry_2d.f90
codee format --config codee-format /tmp/source_registry_2d.f90
diff -u src/numerics/source_registry_2d.f90 /tmp/source_registry_2d.f90
```

Codee is the standard formatter for new or deliberately formatted source, **not
authorization to reformat all existing files at once**. Keep formatting changes
in reviewable groups; avoid unrelated formatting in numerical fixes. Verify a
second format pass makes no further changes and inspect every diff for unchanged
identifiers, literals, physical formulas, floating-point evaluation order,
preprocessor constructs, character continuations, Doxygen header adjacency and
OpenMP directives. Run the existing style/documentation checks, builds and
affected numerical regressions after formatting.

Codee's blank-line rules only control consecutive empty lines and separation
between procedures. They **cannot** identify conceptual sections in scientific
algorithms. Continue to place one blank line between declarations,
initialization, conversion, source evaluation, admissibility checks,
commit/rollback and cleanup when appropriate. Maintain the previously required
blank-line grouping around major control-flow constructs without separating
CASE/ELSE labels from their constructs or Doxygen headers from procedures.
Doxygen coverage, Doxygen line wrapping at 100 columns, lexical correctness
and scientific invariants require dedicated checks and human review.

The formatter is not an automatic CI rewrite step; retain the existing style
and Doxygen checks. Do not make formatting mandatory for every pre-existing
file until representative modules and idempotence have been validated.

## Canonical names and scientific ownership

- Give the same physical quantity the same name through producer, caller and
  callee when meaning, units, representation and time level are identical.
- Use `rho_m`/`inv_rho_m` and `rho_c`/`inv_rho_c` for mixture/carrier density
  and its reciprocal in both REAL and COMPLEX physical states.
- Use `flow_depth`, `total_solid_fraction`, `solid_density`,
  `carrier_density`, `time_step` at mathematical-kernel interfaces when they
  describe those physical quantities. Do not rename unrelated local algebra.
- Use the established global layout vocabulary `n_vars`, `n_eqns`,
  `n_add_gas`, `n_stoch_vars`, `n_pore_vars` in `model_layout_type`.
  Distinguish primitive indices (`idx_h`, `idx_T`, `idx_u`, `idx_v`)
  from conservative equation indices (`idx_totMassEqn`, `idx_uEqn`,
  `idx_vEqn`, `idx_engyEqn`), even when integer values coincide.
- Name each owned configuration by its domain: `rheology`, `thermal`,
  `pore_pressure`, `mass_exchange`. A kernel-local `config` is reasonable
  when the procedure is exclusively about one explicitly typed law.
- Never conflate `q` with `qp`, individual layers, pre/post-step values,
  horizontal/total/squared speeds, or fluid vs mixture density.
- Renames must not change an arithmetic expression, its evaluation order,
  mathematical constants, NAMELIST group names or physical thresholds.
