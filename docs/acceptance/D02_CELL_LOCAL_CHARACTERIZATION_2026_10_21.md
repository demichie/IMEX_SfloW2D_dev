# D02 — Cell-local numerical characterization

Date: 21 October 2026. Baseline: `5aee5d9e` (D01 inventory commit).
Status: **IMPLEMENTED, ACCEPTANCE GATE OPEN**.

The characterization harness and local evidence are complete for the matrix
below. Production defects and configuration/acceptance questions remain.
This is not a successful numerical acceptance and must not be described as
“all D02 tests pass”. No production Fortran, input, RK tableau or tolerance
has been changed. The existing stage test driver gains optional D02 fixtures;
its default behavior is preserved. Standalone checks are included in the
source distribution, without adding a known-failing gate to CI or `make check`.

The [summary JSON](results_d02_cell_local_2026_10_21.json) records provenance,
counts, failures and representative measurements. The [test README](../../TESTS/TEST_RHEOLOGY_CELL_LOCAL/README.md)
documents exact commands, dependencies, norms and tolerances. Full per-case
namelists, logs and JSON results are retained locally under `build-local/d02/`.
They are ignored artifacts, not required files for running the tests afresh.

## Scope and validation design

Classification: numerical characterization/test infrastructure, T1 plus focused
T2 integration. Reference source is the unchanged code audited in D01. There
is no migration or correction of physical formulas in this change set.

| Surface | Coverage |
| --- | --- |
| Laws | Active 1–6, 9–12; disabled rheology; separate legacy-8 ODE probes |
| State/directions | x, y, oblique motion, rest, exact dry and thin states; slope; curvature and excess pressure for laws using them |
| Lava/lahar | Temperature-dependent drag, separate/combined radiation and convection, requested ambient temperature, h threshold, Tc=20 C from both sides |
| Stochastic laws | Frozen Z; lower/upper xi clips; shifted/clipped Froude and h=0.01; Lucas h=0.1 and U=U_w switches |
| Granular laws | Mono/multi-class mu(I), effective load and pressure limits |
| Pore diffusion | P negative/zero/positive; all four inhibition modes; static/dynamic permeability; mono/multi-class solids; Sutherland off/on; packing limits |
| Pore hazards | Rheology disabled with pressure enabled; dry dynamic permeability; zero-solid mono/multi-class dynamic permeability |
| Gas loss | Sixteen isolated Darcy-loss fixtures; independent volume, mixture mass, momentum, energy and transported-pressure balances |
| Newton | Source CS/FD comparisons, production normalized perturbation scaling, eight actual local stage solves and independent residual recomputation |
| Production integration | Quadratic DIRK, Coulomb stopping and coupled lava cooling/drag through the actual integrator; 2 RK stages, 1/4 threads, optimized/strict-debug |
| Actual input parser | Valid control, active model 0, unsupported 7, legacy 8, missing/zero/negative law parameters |

Both source profiles use GNU Fortran 15.2.0 on macOS ARM64. Strict modules use
`-O0 -g -fcheck=all -finit-real=snan`; optimized modules use
`-Ofast -g -funroll-all-loops`. Both test executables enable invalid/zero/overflow
traps. The driver and finite-difference subtraction remain strict. This is a
serial cell-local matrix, not an ownership/thread-safety certification.
The focused production integration separately exercises actual teams of 1/4.

The independent Python oracle constructs conservative states from volume
inventories and evaluates the source equations without consuming production
decoded states. REAL/COMPLEX value parity is checked separately. It preserves
the known baseline ambient overwrite in its formula transcription, while a
separate configuration-fidelity check fails on that overwrite.

Physical-coordinate CS uses relative amplitude `1e-25`; finite differences
sweep `1e-5,1e-6,1e-7,1e-8`. Agreement must hold at two successive step sizes.
The production normalized probe uses `n_vars*epsilon(wp)` and
`q_org=max(abs(q),1e-3)`. Explicit mass/composition columns and currently selected
Newton columns are distinguished. Zero pressure uses its positive-side FD.
Values at named nonsmooth switches are checked against the relevant one-sided
limits instead of asserting a unique derivative. See the README for error
denominator floors, which are fixed in the harness.

## Results

| Check | Result |
| --- | --- |
| Strict cell-local fixtures | **221/249 PASS, 28 FAIL** |
| Optimized + traps cell-local fixtures | **213/249 PASS, 36 FAIL** |
| Production integration, optimized, 1/4 threads | **6/6 PASS** |
| Production integration, strict-debug, 1/4 threads | **6/6 PASS** |
| Actual optimized parser checks | **5/10 PASS, 5 FAIL** |
| Independent gas-loss checks | **16/16 PASS in each profile** |

Counts include isolated path diagnostics; four unique optimized SI rest faults
appear both in combined fixtures and in their isolated SI reproductions.
They are not eight distinct physical defects. Repeated identical counts across
future revisions are not required: the suite has no expected-failure list.

Among completed finite evaluations, the maximum independent source-value error
is `5.06e-12` strict and `6.53e-12` optimized (limit `2e-10`). Unperturbed
REAL/COMPLEX value errors are zero in this matrix. Normalized CS column scaling
errors are below `1e-13` (limit `2e-10`). These value/scale checks do not imply
correct source derivatives.

The source-only integration fixtures check the central cell of a uniform flat
3x3 initial state for one step. It has zero initial spatial residual; no old
Python spatial oracle is involved. Independent roots cover quadratic DIRK and
thermal balance; Coulomb stopping is checked exactly. Known stage assembly,
the stored effective increment after corrections, raw final assembly and final
projection agree in the tested cases. The largest source-stage state error is
`2.90e-6` for quadratic drag, within the configured solve accuracy `1e-5`.
The coupled lava error is `5.61e-11`; Coulomb stopping is exact. Central states
are bitwise identical between 1 and 4 threads for each tested build.

## Findings requiring resolution before acceptance

### D02-01 — Smooth source derivatives are incomplete

Twenty fixtures in each profile fail CS/FD agreement:

- Sixteen STATIC/DYNAMIC inhibition cases lose dependence through the REAL
  polynomial argument and/or REAL height-dependent transition. The discrepancy
  includes conservative energy column 4, selected by Newton in the model-3
  coupled-temperature diagnostic.
- Four multi-class dynamic-permeability cases with OFF/STEP inhibition fail
  composition columns 5/6. Their Sauter diameter/sphericity calculations use
  REAL volume fractions. Those composition columns are explicit in the current
  Newton mask; this distinction is retained in the results.

Representative strict case `pore-2-DYNAMIC-0-1`, pressure-source row versus
conservative energy column: CS `+5.90186871499351e-4`, finest FD
`-1.29728696917244e-3`. The sign is wrong and the plateau discrepancy is about
`1.455`; this is not finite-difference cancellation or a switch derivative.

The coupled-temperature probe uses gas-carried material with selector 3 to
exercise the API's thermal/pressure coupling, a combination the input selector
does not prohibit. It is a synthetic derivative fixture, not a claim of a
physically calibrated gas-lava application.

**Closure evidence:** retain complex dependence through smooth intermediate
algebra (or explicitly define/freeze a different supported mathematical
contract), preserve REAL source values, and rerun the regular/one-sided matrix
plus relevant local Newton/integration tests. Do not hide the issue by removing
composition checks or changing tolerances.

### D02-02 — Initialization and singular diffusion coefficients

Four isolated hazard cases fail in both profiles:

- `pore-only`: pressure is enabled with rheology disabled; the source skips
  `c_phys_var` but subsequently uses its local state outputs.
- `dry-dynamic-pore`: coefficient evaluation occurs before the wet-state guard.
- `zero-solids-dynamic-pore` and `zero-solids-mono-pore`: dynamic permeability
  divides by vanished solid inventory; multi-class weights also have zero sum.

The strict profile initializes unassigned local REAL storage to signaling NaN,
so the first case cannot accidentally pass by inheriting finite stack values.
The optimized profile also traps on these recorded inputs. This supersedes the
static-only status of the corresponding D01 findings.

**Closure evidence:** define supported flag/phase combinations; initialize
physical state independently of the rheology flag where pressure requires it;
guard dry/zero-solid coefficient evaluation or reject unsupported configurations
before entering the kernel. Test the chosen limits without manufacturing
physical carrier loss or an arbitrary permeability floor.

### D02-03 — Input validation admits missing/zero denominators

The actual parser admits five configurations this validation check requires it
to reject: model 1 with missing xi; model 9 with missing mu_0 or zero Fr_0;
model 10 with zero U_w or missing mu_inf. The short solver runs reach their end,
so an unrelated startup error is not being confused with rejection.

The valid control and legacy 8 are accepted; active 0, unsupported 7 and the
all-negative model-9 input are rejected as expected. Rejections use Fortran
STOP with exit code zero on this toolchain. The checker therefore asserts the
specific rejection diagnostic, not merely a nonzero process status.

**Closure evidence:** validate each required parameter individually, distinguish
permitted zero-friction coefficients from strictly positive denominators and
retain valid endpoint cases. Changing selector names/removing 8 remains D08/D04.

### D02-04 — Optimized SI rest arithmetic with traps

Models 5, 6, 11 and 12 trap in combined rest fixtures only in the optimized
profile. Isolated reproductions localize the fault to the REAL SI path; the
REAL implicit and zero-perturbation COMPLEX implicit paths return zero.
Strict SI calls return zero too.

For these four inputs, the identical optimized production objects linked to a
test main without FPE traps return finite zero SI vectors. Thus the observation
is **optimized arithmetic exposed by trapping**, not demonstrated failure of a
normal untrapped production run. `-Ofast` permits transformations that do not
preserve trapping behavior; compiler flags and guarded zero-speed expressions
both require review. No successful untrapped call replaces a trapping failure
in the acceptance counts.

**Closure evidence:** isolate the generated operation/flag interaction and adopt
a consistent supported build/trapping contract; verify rest, dry and moving
states under that contract, including another toolchain where relevant.

### D02-05 — Configurable ambient temperature is overwritten

Requesting `T_env=350` yields `T_env_after=300`. The REAL source matches the
historical hard-coded cooling formula, but input configuration does not survive
the source evaluation. This is a shared write as well as a thermal policy issue.

**Closure evidence:** decide whether the ambient must be configurable or fixed,
make initialization/input/kernel behavior consistent, and verify cooling values
and derivatives for ambient/reference/flow temperatures. Ownership belongs in
the D03/D06 migration; D02 must not silently certify the current override.

### D02-06 — Convergence does not guarantee the residual threshold

The local solves for lahar, quadratic drag and pore diffusion return converged
while the independently recomputed stage residual exceeds
`tol_abs+tol_rel*abs(initial_residual)`. Maximum component ratios to this
threshold are approximately `9.03`, `1.35`, and `1.32`, respectively.

The current solver also permits convergence on a small normalized Newton
correction before applying that correction. This finding does not claim the
existing code promises residual-only convergence; it shows that its alternative
exit criterion does not guarantee the audited residual threshold. The independent
quadratic stage state remains within the current relative solve accuracy.

**Closure evidence:** document and accept a justified stopping contract or
require an updated state to satisfy the chosen residual criterion. Preserve
residual/error measurements and revalidate scaled Newton and coupled sources;
do not simply increase the test's residual bound to obtain PASS.

## Baseline differences retained for future extraction

Both pressure-diffusion and Darcy-loss values match their respective independent
formulas. They use different multi-class sphericity weights (`alpha*d^2` versus
`alpha/d`) and different dry/solids guards. Their successful individual tests
do not establish that the two coefficients can be merged into one provider.
This remains a scientific/ownership decision for D06 with explicit comparisons.

The legacy model-8 moving/rest/dry ODE probes match their independent speed-decay
formula. That does not turn the method into a DIRK stage solve. D04 removal is
unchanged. Model 7 remains unsupported. Global xi scratch and same-value writes
at common parameters are not certified thread-safe by the homogeneous T2 runs;
D03 still requires an ownership/access audit.

## Remaining work and gate boundary

To close D02, resolve D02-01/02 and input-domain failures, explicitly settle the
trapping/thermal/stopping contracts in D02-04/05/06, and rerun the corresponding
targeted T1/T2 checks. Record independent source value/derivative equivalence
when changing algebra. Cases cannot be waived solely because previous broad
integration tests passed.

There is no evidence here for all production reologies in real terrain, every
mass-exchange/deposition configuration, all prior Python suites, Linux/Docker,
Intel ifx or MPI. These were not rerun: production formulas are unchanged and
this characterization targets the local kernels and a small actual-stage
integration matrix. No GitHub Action, commit, push or release was performed.
