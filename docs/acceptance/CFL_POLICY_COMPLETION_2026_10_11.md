# CFL policy implementation and acceptance

The authorized timestep contract offers `directional_min` as the unchanged
default and `combined` as an explicit option. The scientific baseline N9
requires a separate approval; this implementation does not change IMEX tableaux,
admissible positive CFL limits, event clipping, growth limiting or retry policy.

## Spatial proposal and accepted step

Both options use the same final reconstructed HP-PCCU face-speed maxima.
`directional_min` retains the previous multiplication/division order.
`combined` uses `CFL/(lambda_x/dx + lambda_y/dy)` in 2D; collapsed x/y grids
use the corresponding legacy 1D expression to preserve exact equivalence.
Zero propagation returns `MAX_DT` without division by zero.

Input accepts precisely `CFL=-1` in addition to the existing positive interval.
This sentinel returns the spatial cap without reconstructing face speeds.
Unknown policy names and other negative CFL values fail explicitly. The
selected policy is stored in effective input backups and invocation manifests.

The integrator still caps the spatial proposal at output/runout/probe events,
the final time and the previous two-step growth bound. These formulas are not a
general stability theorem for the supported `N_RK=2,3,4` tableaux. Current-IMEX
common-time comparisons and historical SSPRK2 comparisons keep their distinct
meaning and unchanged limits; universal autonomous roundoff agreement is not
asserted by selecting the same spatial CFL formula.

## Direct policy and temporal verification

The spatial fixture checks both policies against one common face-speed payload
and the pinned independent Python combined bound. Rotated x/y grids, unequal
spacings, a one-cell domain, an empty workset, caps and the sentinel are checked
with strict/optimized production objects and actual one/four-thread teams.

The production temporal fixture covers all three stage counts, omitted/default/
combined selections and the sentinel. Independent time algebra checks require
`64*eps64` agreement and non-vacuous growth, three output events and a final
event. A test-only observational main must preserve every conservative output
and binary checkpoint byte for byte.

The preceding `691a3e21` executables are separate historical comparators.
For each profile/team/stage combination, uninstrumented fields/checkpoints and
independently observed spatial proposals/accepted steps must match exactly
when the new parameter is omitted. The historical observer must also preserve
its own canonical outputs.

## Retained diagnostics

The first collapsed-axis check exposed a one-bit difference between equivalent
division orders, fixed only in the new combined 1D branch. The default branch
was not changed.

The first temporal fixture reused the compact, open-boundary Gate-H domain.
At `t=0.095 s`, mass reached its north/east/south boundary and the apparent
closed-domain budget failed: mass-sum drift `-5.809124559164047e-6`, bound
`1.0913936421275139e-6`. This is retained as a fixture diagnosis, not waived.
Only the new temporal fixture adds eight dry cells per side; the same budget
bound is retained. Frozen Gate-H generation defaults to zero padding.

The obsolete first-candidate audit was interrupted after retaining its failure
logs. Its partial results are not counted as a completed acceptance audit.

## Validation status

The [pinned completion evidence](results_cfl_policy_2026_10_11.json) records:

- **45/45 full-audit groups pass**, with eight passing cross-thread pairs.
- All 16 preceding full-solver groups retain bitwise-identical canonical
  conservative fields and checkpoints under the default policy.
- 128 independent spatial comparisons and 12 collapsed x/y/point runs pass.
- 72 temporal cases cover all policies/sentinels/stages in both profiles and
  actual one/four-thread teams. All 36 temporal thread pairs retain identical
  fields, checkpoints and proposed/accepted steps.
- 12 historical stage/profile/team comparisons retain bitwise-identical
  proposed/accepted times; their observational mains preserve canonical outputs.
- Eight additional original-duration `combined` deterministic/stochastic
  restart groups and all four corresponding bitwise thread pairs pass.

The production-source hashes match the final audit snapshot. The subsequently
added historical-time helper and explicit capture-adapter policy override are
separately verified and pinned; neither changes numerical sources. The closure
metadata/checker updates are verified separately after numerical completion.
Source documentation and the CFL, Gate-H and geometry checker controls pass.

D-N7-CFL and **N7 are satisfied within these defined contracts**; N8 remains
satisfied. **N9 remains open**, with no scientific-baseline approval implied.
The old autonomous Python/Fortran diagnostic is retained as not satisfied:
choosing the same spatial formula is not a claim of universal equality between
independently adapted IMEX trajectories. Common-time field/stage gates and
historical SSPRK2 error limits retain their separate interpretation.
