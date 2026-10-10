# Dynamic geometry and excavation field diagnostics

The diagnostic history below preserves the failed finite-time geometry and
all-cell excavation comparisons and the mechanisms isolated before each
authorized correction. The final corrected candidate now passes 41/41 audit
groups and closes N7-C/N7-D in their defined fixtures; see the
[completion report](N7_CD_COMPLETION_2026_10_10.md). The separate CFL decision
and N9 approval remain open. None of the earlier discrepancies was accepted
by changing tolerances.

## Comparisons and retained limits

The current-policy comparison uses common IMEX tableaux, initial conservative
fields, Q1 bed, ambient density and time-step policy. Its all-cell thickness
and free-surface limit is
`32768 * float64_epsilon * accepted_steps * max(1, initial_h_max)`.
For the original coarse Gate-H trajectory this was approximately
`9.459e-10 m`; the bound scales with the actual accepted step count.

The separately regenerated historical SSPRK2 trajectory keeps its own
harmonic CFL and stage-retry policy. The user approved normalized L1 limits
of `2e-4` for thickness and `1e-4` for free surface and a maximum pointwise
difference of `0.003 m`. Existing aggregate and uphill limits are unchanged.
Old scalar port-audit norms are not substitutes for actual reference fields.

G=1 fields passed all three grids in the focused strict run. The maximum
thickness differences were approximately `6.153e-12`, `1.245e-11` and
`3.422e-11 m` at spacings `0.4`, `0.2` and `0.1 m`. Historical SSPRK2
differences also passed the approved limits in these cases.

For slope plus curvature at spacing `0.4 m`, explicit matching of production
LS boundary-strip copying reduced the discrepancy from approximately
`2.985e-4 m` to `5.272e-5 m`. Closing final external ghosts and applying the
existing exact-rest pressure roundoff rule did not remove the remaining
trajectory discrepancy. These are reference-policy distinctions, not
changes to the checksum-pinned clean kernels or relaxed acceptance limits.

## Isolated amplification mechanism

A test-only trace build recorded conservative states before and after all
65 accepted IMEX steps. Its final canonical output is byte-identical to the
unmodified solver output. Starting the independent reference afresh from
each actual production state gave maximum one-step differences of
`2.274e-13 kg/m2` in mass and `5.116e-13` and `5.143e-13` in the two momenta.
Thus the large accumulated field discrepancy is not a comparably large
one-step spatial or temporal assembly error in this trajectory.

The first different exact-rest face decisions occur before step 21, at
`t=0.00427287130950639 s`, on x faces 24 of zero-based rows 4 and 17.
Before that step the maximum thickness difference is only
`3.865e-15 m`. The nominally dry neighbouring cell is exactly zero in
production but retains approximately `1.712e-17 m` in the independent
trajectory. Its desingularized x velocity is approximately `0.002581 m/s`.
Although both depths are below the existing face dry tolerance, the raw
auxiliary centre velocity participates in direct reconstruction before the
final safe-velocity momentum bounds.

The wet endpoint on the first differing face has depth approximately
`4.175e-5 m` in both trajectories. Its reconstructed volumetric momentum is
exactly zero in production and approximately `9.245e-8 m2/s` in the
reference. The exact-rest scalar guard therefore blocks mass transport in
one trajectory but not the other. The reference unfiltered face mass flux
is approximately `2.844e-4 kg/(m s)`. The next state develops a maximum
thickness difference of `3.229e-10 m`; the difference then accumulates.

This identifies a sensitivity of the reconstruction and exact-rest transport
switch to auxiliary velocities in nominally dry neighbouring cells. A
candidate correction is to use zero auxiliary velocities in reconstruction
at the existing dry depth threshold while retaining conservative cell mass.
That would be a numerical policy change requiring its own before-and-after
validation, not merely an optimization or a new comparator tolerance.

## Reproducibility and status

The versioned dynamic contracts preserve the original failed quiet-collar
and circular-refinement assessments. The user-approved circular measure
separates initial sampling from the fourfold moment introduced by evolution,
without changing the absolute final bound.

Focused failure evidence is retained in
`/tmp/imex-pccu-gate-h.Naj6sp` and the trace diagnosis in
`/tmp/imex-gate-full-trace.Hh0s8u`. Machine-readable diagnostic summaries
are recorded in the accompanying JSON. Temporary paths are local evidence,
not portable checked-in golden fields.

The production solver was unchanged throughout the trace diagnosis above.
The user subsequently authorized the dry auxiliary reconstruction correction.
Its first candidate passed the 69 dynamic cases in both compiler profiles,
but did not close Gate-H: the coarse slope-plus-curvature field difference
was approximately `1.159e-5 m`, above the unchanged `2.299e-9 m` limit, and
the next refinement stopped on floating-point overflow. The debugger located
the overflow in the squared velocity used by physical-state conversion.

The local curvature force still used unresolved centre velocities below the
shared dry depth threshold. The user therefore also authorized a guard on
that explicit source only: evaluate curvature for
`h > dry_thickness_tolerance`, retaining positive conservative mass, other
sources, the existing CFL and all acceptance limits. Direct caller tests
cover depths below, exactly on and above the threshold; the current-policy
reference specifies this guard separately from the unchanged historical
SSPRK2 core. `contract_v4.json` preserves the first candidate's scope and
`contract_v5.json` records the additionally authorized source guard.

With both corrections, all strict focused grids complete without the
previous overflow. Independently adapting accepted times still fails the
roundoff comparison: thickness differences are approximately `6.747e-9`,
`3.256e-5` and `6.253e-9 m` on the three refinements. These are retained
failures, not accepted by loosening the limit.

The observational coarse trace has 159 steps and produces byte-identical
canonical output. The two autonomous time sequences differ by up to
`1.644e-6 s` in accepted dt. Replaying the independent trajectory from the
common initial fields at the actual accepted times reduces the final
thickness difference to `2.240e-11 m`, below the unchanged `2.314e-9 m`
bound. No actual intermediate state is used as a reference field in this
common-time evolution.

The user explicitly approved separating the common-time operator/stage gate
from the autonomous adaptive-control comparison. `field_contract_v2.json`
records that separation, without changing either field limit. The latter
comparison remains non-satisfied under `D-N7-CFL`; it does not become a
passing assertion by omission. The complete common-time audit with the
auxiliary and curvature corrections finished with 37 of 41 groups passing.
The intermediate slope-plus-curvature grid still failed in both profiles
and thread teams: approximately `3.256e-5 m` against the unchanged
`2.678e-9 m` strict bound. All other current-policy grid comparisons and all
historical field comparisons passed. This did not close N7-D, and N9 remains
open.

## Dry cell momentum memory on the intermediate grid

In the strict intermediate trajectory, zero-based cell `(row 37, column 46)`
retained depth `2.196349e-15 m` from steps 6 through 47. Its conservative
momenta nevertheless grew to approximately `5.110335e-6` and
`8.134226e-6 kg/(m s)`. Auxiliary reconstruction velocities were already
zero and the local curvature source was already disabled in this dry cell.
The independent trajectory had an effectively empty cell and did not retain
those moments.

Upon rewetting, the cell reached depth `2.021297e-7 m`. Its retained moments
then affected resolved velocity bounds and neighbouring wet-cell fluxes.
The maximum thickness difference grew from `5.777e-12 m` before step 48 to
`1.315e-8 m` before step 49 and exceeded the roundoff gate. At step 48 the
normal exact-rest face masks agreed in both directions, and the curvature
depth masks agreed throughout the trajectory. This diagnosis is distinct
from the earlier coarse-grid exact-rest face-switch discrepancy.

Fresh single-step replays from sampled actual states agreed within
`2.274e-13 kg/m2` in mass and `7.854e-13 kg/(m s)` in momentum. Such local
replays diagnose the amplification mechanism; they do not replace the
independent finite-time field gate. The trace canonical output was
byte-identical to the uninstrumented output. Trace and replay evidence are
retained in `/tmp/imex-guard-common-dx0p2.5st_mqop`.

The user authorized a further conservative-state correction: reset only
the two moments at `h <= dry_thickness_tolerance`. The implementation
applies it to the accepted final state after raw observation and checks,
retains mass, thermal energy and every transported component, and updates
the corresponding primitive momenta and auxiliary velocities. It does not
change raw IMEX stages, implicit solves, CFL, thresholds or limits.
`contract_v6.json` and `field_contract_v3.json` preserve this added scope.
The new complete audit passes 41/41 groups and all six cross-thread
Gate-H/restart comparisons. The corrected intermediate common-time thickness
error is `6.180e-12 m`, below the unchanged `2.183e-9 m` bound. The new trace
retains positive dry mass while all accepted dry momenta are exactly zero.
The supplementary autonomous comparison also passes on this strict
intermediate case (`4.052e-10 m`), but it is not a hierarchy-wide approval of
the adaptive-control contract. N7-C and N7-D close; `D-N7-CFL` remains open.
