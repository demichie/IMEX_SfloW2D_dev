# Fortran source style (mandatory)

Applies to all production `src/**/*.f90` and to Fortran tests when edited.

## Indentation and layout

- Four ASCII spaces per lexical nesting level; **no TAB characters**, including in comments.
- Align `END IF`, `END DO`, `END SELECT`, `END BLOCK` and `END SUBROUTINE name` with their respective openers.
- Indent continued statements consistently with spaces; preserve the semantics of free-form `&` continuation and OpenMP/compiler directives.
- Separate declaration blocks, initialization, state conversion, source evaluation, admissibility checks, commit/rollback and cleanup with blank lines. Avoid stacked unrelated statements.
- Separate procedures visibly with a blank line and their own Doxygen header; use meaningful section comments rather than decorative banners inside procedures.
- No trailing whitespace. Pure formatting changes must preserve all non-whitespace source tokens, floating-point expression order, runtime behavior and OpenMP directives.

## Mandatory Doxygen coverage

Document every module, derived type and significant field or variable (meaning, units, rank, ownership when relevant). Every FUNCTION and SUBROUTINE, including PRIVATE, PURE, ELEMENTAL, internal and type-bound procedures, must have an immediately preceding `!> @brief` header, `@details` when needed, and **one correctly directed** `@param[in]`, `@param[out]` or `@param[in,out]` for every dummy argument, including `this`; functions also require `@return`. Describe dimensions, units, preconditions, error behavior, side effects and thread-safety where applicable. Never invent generic documentation to satisfy a checker.

Comment substantial algorithmic blocks with the reason for ordering and key invariants. Do not narrate trivial assignments.

## Checks

The repository already provides `python3 TESTS/TEST_SOURCE_DOCUMENTATION/check_documentation.py` for Doxygen contract checks, including argument names/directions; it also supports `--check-against-git REF` to compare executable lines. The complementary `python3 tools/docs/check_fortran_style.py` checks four-space indentation, TABs and trailing whitespace. During migration, run it with explicit paths for the files normalized so far. Run it over `src/` only when all 64 production sources have been normalized.

This is a formatting/documentation phase, not an opportunity to alter formulas or the Lie algorithm.
