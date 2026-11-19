# Mu(I) analytic quadrature failure diagnostics

This is an opt-in diagnostic for the existing single-layer `EXAMPLE_DAMBREAK`
and FINH cases. It does not change any physics, numerical methods, tolerances,
integration limits or reconstruction.

Generate the dam-break input as usual and set
`RHEOLOGY%RHEOLOGY_METHOD = 'ANALYTIC'` inside `&RHEOLOGY_PARAMETERS`.
Run from the example working directory with one OpenMP thread:

```sh
IMEX_LIE_DIAGNOSTICS=1 OMP_NUM_THREADS=1 ./IMEX_SfloW2D > solver.log 2>&1
```

Use the executable path appropriate for your local setup. The error/witness
records are written to standard error. The final `MUI_ANALYTIC_CELL` record
preceding the failure gives the cell indices, layer, Lie group, source time and
source timestep; `MUI_ANALYTIC_Q` is the conservative input vector.
`MUI_ANALYTIC_STATE` records depth, mixture density, horizontal velocity,
timestep, speed ratio and reduced gravity. `MUI_ANALYTIC_STRESS` records
effective normal stress, curvature acceleration and initial speed.
`MUI_ANALYTIC_COEFFICIENTS` lists the five frozen coefficients supplied
to the rational mu(I) exact-flow law in their original order.

At an adaptive Simpson failure the solver also emits
`EXACT_QUADRATURE_FAILURE` (flow model, interval, panel estimates and delta),
`EXACT_QUADRATURE_TOLERANCE` (requested tolerance, observed absolute delta,
and acceptance threshold) and `EXACT_QUADRATURE_COEFFICIENTS`.
These failure-only records are emitted even if Lie diagnostics are disabled.

With multiple OpenMP threads cell witnesses can interleave; use one thread
for attribution to a single failing cell. Do not change the quadrature depth
or tolerance to make the run pass. Attach the final witness and failure
records for numerical investigation.
