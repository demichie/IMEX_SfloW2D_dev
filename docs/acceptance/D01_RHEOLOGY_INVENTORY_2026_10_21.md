# D01 — Rheology and cell-local source inventory

Date: 21 October 2026. Audited revision:
`65888eeceacbf78497bb3f7d361dc5e598d8ef2d`.

Status: **static inventory complete**. This report describes the implementation,
not an endorsement of its physical formulas or a numerical acceptance of D02.
No source, input, RK tableau or numerical tolerance is changed by D01.
The companion [JSON inventory](results_d01_inventory_2026_10_21.json) records
source hashes, routine spans, dispatch coverage, parameters and the D02 backlog.

## Scope and evidence

The audit follows the input selector through source evaluation, state decoding,
Newton differentiation, IMEX stage corrections, mass exchange and diagnostics.
Names attributed to laws below are the names used by the code/roadmap; no
literature-to-code equivalence is asserted. The algebra below is a transcription
of the code, including its guards and historical conventions.

Primary entry points (paths are relative to the repository root):

| Responsibility | Current implementation |
| --- | --- |
| Select/validate laws; read and derive parameters | `src/io/inpout_2d.f90`: `init_param`, `read_param`; namelists `rheology_parameters`, `temperature_parameters`, `pore_pressure_parameters` |
| Declare material/rheology coefficients | `src/core/constitutive_parameters_2d.f90` |
| Declare flags, selector and indices | `src/core/parameters_2d.f90` |
| Assign implicit equation mask | `src/physics/equation_terms_2d.f90`: `init_problem_param` |
| Evaluate implicit sources, REAL and COMPLEX API | Same module: `eval_implicit_terms` |
| Evaluate semi-implicit friction, REAL only | Same module: `eval_nh_semi_impl_terms` |
| Special model-8 ODE integrator, REAL only | Same module: `integrate_friction_term` |
| Recover physical state and mixture properties | `src/state/state_conversion_2d.f90` and its shared `.inc` kernels |
| Evaluate stage residual/Jacobian | `src/numerics/nonlinear_solver_2d.f90`: `eval_f`, `eval_jacobian`, `solve_rk_step` |
| Apply SI/Newton corrections and store effective RK terms | `src/numerics/time_integration_2d.f90`: `imex_RK_solver` |
| Apply deposition, carrier loss and thermal mass exchange | `src/physics/equation_terms_2d.f90`: `eval_mass_exchange_terms`; `src/physics/mass_exchange_2d.f90` |
| Supply stochastic forcing | `src/stochastic/stochastic_module.f90` |
| Write law-dependent diagnostics | `src/io/inpout_2d.f90`: `write_netcdf_timestep` |

These are current owners/responsibilities, not the symbol-by-symbol ownership
manifest required by D03. Future `rheology_*.f90`, `thermal_laws_2d.f90` and
`pore_pressure_2d.f90` do not exist yet.

## State and notation

For the local conservative state, `M=q(1)=rho*h`, `p=q(2:3)`, and
`E=q(4)=M*cp_mix*T`; E contains thermal energy, not kinetic energy.
Solid and additional-gas conservative entries are component masses per area.
The transported pore variable is `q(idx_pore)=M*P`, with P the excess pressure.
The implicit source writes `idx_poreEqn`; the primitive entry holds P itself.

Let `v=(u,v)`, `Vh=sqrt(u^2+v^2)`, and
`V=sqrt(u^2+v^2+w^2)`. The code uses `w=u*Bx+v*By` when slope correction is on,
otherwise zero. `g'=g*(rho-rho_a_amb)/rho`, `G=grav_coeff`, and
`C=u^2*Bxx+2*u*v*Bxy+v^2*Byy` is the signed curvature acceleration.
`A=sum(alpha_s)` denotes total solid volume fraction, not mass fraction.
Negative source vectors below are added to the conservative evolution.

The implicit API accepts either a paired REAL input/output or a paired COMPLEX
input/output. Both execute the **same COMPLEX source algebra**: the REAL input
is promoted to COMPLEX with zero imaginary part, and its result is projected
back to REAL. This is not a separate REAL implementation of each implicit law.
The semi-implicit evaluator and model-8 integrator have no COMPLEX API.

## Selector and model inventory

With `RHEOLOGY_FLAG=T`, `read_param` accepts exactly
`1,2,3,4,5,6,8,9,10,11,12`; zero, 7 and unknown IDs reach an error/STOP branch.
With the flag off, zero is the conventional no-rheology selector, not an
accepted active law. Do not infer support from the range 0–12.

| ID | Implemented law / planned name | Implicit source | Semi-implicit source | Law parameters (besides common state/material data) |
| --- | --- | --- | --- | --- |
| 0 | Disabled / `none` | None when flag is off | None when flag is off | `rheology_flag` |
| 1 | Voellmy–Salm / `voellmy` | Turbulent drag | Coulomb load | `mu`, `xi`; stochastic Z |
| 2 | Plastic / `plastic` | Constant-strength drag | Empty branch | `tau` |
| 3 | Temperature-dependent lava / `lava` | Viscous drag and radiation/convection | Yield strength | `nu_ref`, `visc_par`, `T_ref`, `tau0`; thermal coefficients |
| 4 | O'Brien/FLO2D lahar / `lahar` | Viscous and turbulent-dispersive drag | Yield stress | `alpha1_ref` (input-local), `alpha1_coeff` (derived), `beta1`, `alpha2`, `beta2`, `Kappa`, `n_td`, `n_td2`, `T_ref` |
| 5 | Kurganov–Petrova Example 5 / `empirical` | Depth-dependent empirical drag | No branch | Hard-coded `1e-3`, `10` |
| 6 | Bursik–Woods quadratic drag / `quadratic` | Quadratic drag | No branch | `friction_factor` |
| 7 | **Unsupported** | No branch | No branch | No validated parameter set or future file name |
| 8 | Special quadratic ODE path / removal in D04 | Bypasses Newton; no source branch | No branch | `friction_factor` |
| 9 | Froude-dependent Coulomb / `froude` | Empty branch | mu(Fr) | `mu_0`, `mu_inf`, `Fr_0`; stochastic Z |
| 10 | Lucas velocity weakening / `lucas` | No branch | mu(U) | `mu_0`, `mu_inf`, `U_w`; stochastic Z |
| 11 | Regularized mu(I) / `mui` | No branch | Granular friction | `mu_s`, `mu_2`, `muI_inf`, `I_0`, `diam_s`, `rho_s` |
| 12 | mu(I) + Voellmy / `mui_voellmy` | Same turbulent drag as 1 | Same mu(I) block as 11 | Parameters of 11 plus `xi`; stochastic Z affects xi |

An empty/no rheology branch does **not** imply the entire implicit vector is
zero if pore-pressure evolution is enabled independently.

### Formulas and guards by law

**1 — Voellmy.** For `REAL(V) != 0`, implicit momentum is
`S_I=-rho*g'/xi_eff * V^2 * v/V`. Stochastic xi is
`xi_eff=clamp(xi+Z,1,1e5)`; otherwise `xi_eff=xi` without that clamp.
`xi_temp` is a shared module variable written during evaluation.
The REAL SI block returns zero for `h<=epsilon(wp)` or `V=0`.
It forms `L=max(0,mu*rho*(g'+C)*h*sqrt(G))` (C only if enabled),
subtracts `mu*sqrt(G)*P` if pore pressure is enabled, clamps again, and applies
`S_SI=-max(0,L)*v/V`. The order of these clamps is part of the baseline.

**2 — Plastic.** For `REAL(V) != 0`, `S_I=-rho*tau*v/V`.
There is no independent SI yield contribution. Preserve the density factor
as implemented; do not silently identify tau with the lava yield parameter.

**3 — Lava.** `gamma=3*nu_ref/max_branch(h,1e-10)*exp(-visc_par*(T-T_ref))`,
where the branch is selected using `REAL(h)>1e-10`. Moving-state momentum is
`S_I=-rho*gamma*v`. Energy is
`S_I,E=-radiative_term_coeff*(T^4-T_env^4)-convective_term_coeff*(T-T_env)`.
The evaluator **assigns shared `T_env=300`** on every model-3 call, overriding
the configurable value. `read_param` derives the coefficients as
`SBconst*emissivity*exp_area_fract` and `atm_heat_transf_coeff*exp_area_fract`.
SI momentum for a moving wet state is `-tau0*v/V`.
`init_problem_param` marks equation 4 implicit for selector 3.

**4 — Lahar.** `Tc=T-273.15`. For `REAL(Tc)<20`,
`alpha1=alpha1_coeff*1e-3*10^a`, with
`a=1301/(998.333+8.1855*(Tc-20)+0.00585*(Tc-20)^2)-1.30223`.
Otherwise `alpha1=alpha1_coeff*1.002e-3*10^b`, with
`b=(1.3272*(20-Tc)-0.001053*(Tc-20)^2)/(Tc+105)`.
Then `eta_fluid=alpha1*exp(beta1*A)`.
`read_param` derives `alpha1_coeff` from `alpha1_ref` at `T_ref` and stores
`n_td2=n_td^2`; active input requires liquid without gas.

For h above the threshold,
`s_v=Kappa*eta_fluid*V/(8*rho*g*h^2)` and
`s_td=n_td2*V^2*h^(-4/3)`.
The threshold is `1e-10` for `wp==sp`, otherwise `1e-20`;
below it the denominators use that threshold.
For `REAL(V)>0`, `S_I=-g*rho*h*(s_v+s_td)*v/Vh`.
SI uses `tau_y=alpha2*(exp(beta2*A)-1)`,
`s_y=tau_y/(g*rho*max_branch(h,1e-20))`,
and `S_SI=-g*rho*h*s_y*v/Vh` when moving.
The same viscosity-temperature formula is repeated in mass exchange, where
it changes settling viscosity and selects the lahar-specific deposition path.

**5 — Empirical.** `c_tau=1e-3*V/(1+10*h)` and, for `REAL(V)!=0`,
`S_I=-rho*c_tau*v/Vh`. The constants have no dedicated namelist parameters.

**6 — Quadratic.** For `REAL(V)!=0`,
`S_I=-rho*friction_factor*V^2*v/Vh`. It uses the ordinary Newton/DIRK path.
Slope-dependent V and the horizontal projection Vh must remain distinct.

**7 — Resolution of the roadmap question.** No active source branch, special
solver path, input validation branch or ID-7 fixture was found. `read_param`
rejects it through the final ELSE when rheology is enabled. The code does not
establish a physical identity for a former model 7: keep it unsupported and
assign no invented law/module name. This finding does not introduce a new
input migration policy or a new rejection diagnostic.

**8 — Special path.** For nonzero mass and `Vh>epsilon(wp)`, the REAL integrator
sets `Vh_new=1/(1/Vh+friction_factor*dt/h)`, scales u/v by `Vh_new/Vh`, then
converts back to conservative variables. `solve_rk_step` calls it with
`dt_step`, declares convergence and returns before Newton, without using
`a_diag` in this update. The integrator later stores the effective increment
divided by `dt*a_diag`. This is an **ODE integration**, not the solution of
the algebraic DIRK stage. The input parser groups 6 and 8 under Bursik–Woods;
NetCDF shear/Rouse diagnostics also have a model-8-only branch.
No ID-8 input assignment was found in the audited tracked test/example files.

**9 — Froude.** The REAL SI block requires wet, moving flow and `g'>0`.
`Fr=V/sqrt(g'*h)`; with stochastic forcing add Z only for `h>0.01`, then
clamp Fr to zero from below. The actual formula is
`muF=mu_inf+(mu_0-mu_inf)*exp(-Fr/Fr_0)`:
its zero/infinite-Fr limits are mu_0/mu_inf, respectively, despite reversed
comments above the branch. `S_SI=-rho*muF*h*G*(g'+C)*v/Vh`;
C remains zero in this block because the routine initializes it to zero and
only other model branches recompute it. Neutral/buoyant states return zero.

**10 — Lucas.** The REAL SI block requires wet, moving flow.
`U_f=V+Z` if stochastic, otherwise V. For `U_f<U_w`, or for `h<0.1` in the
other branch, `muF=mu_0`; otherwise
`muF=mu_inf+(mu_0-mu_inf)/(U_f/U_w)`.
The momentum formula is the same as 9 with this muF, but there is no `g'>0`
guard and C is again zero. Preserve the branch order for negative shifted U.

**11/12 — mu(I).** For moving wet flow, select Sauter diameter and mean solid
density for multiple classes, otherwise use the first class directly.
The code floors h to epsilon, computes C only when enabled, and forms
`sigma_v=max(0,rho*(g'+C)*h-P)` (omit P if disabled).
`sigma_n=max(0,sigma_v*sqrt(G))`, then floors it to epsilon.
`shear_rate=2.5*V/(h*sqrt(G))`,
`I=diam_characteristic*shear_rate/sqrt(sigma_n/rho_particle)`,
`mu_I=(mu_s*I_0+mu_2*I+muI_inf*I^2)/(I_0+I)`.
SI momentum is `-mu_I*sigma_n*v/V`.
The 12 block duplicates 11 and adds the implicit turbulent term of 1.
Both SI blocks are REAL only; neither is differentiated by the current
complex-step Newton Jacobian. A future extraction must not accidentally
reclassify these terms as implicit.

## Thermal and pore-pressure source inventory

### Thermal responsibilities

Only model 3 has a dedicated radiative/convective energy relaxation source
in `eval_implicit_terms`. `thermal_conductivity`, `T_ground`, `c_p`, `emme`
and `enne` are declared/read but have no executable uses in current source
evaluation; their presence in the namelist does not establish conduction,
viscous heating or another active cooling law.

The explicit source evaluator sets the local hydrostatic/curvature energy
component to zero, then adds energy carried by configured external sources.
Mass exchange separately adds/removes thermal energy with deposition,
erosion, carrier losses, entrainment and vaporization. These are not part of
the model-3 cooling ODE and must not be dropped in D06. Its lahar viscosity
duplication must also be tracked across extraction.

### Excess pressure diffusion

With `pore_pressure_flag`, `eval_implicit_terms` computes
`gamma_gas=cp_air/(cp_air-R_air)`, `beta_gas=1/(gamma_gas*pres)`,
`rho_gas=pres/(R_air*T)`, `porosity=1-A`, and
`D=k/(porosity*nu_carrier*rho_gas*beta_gas)`.
For `REAL(P)>=0` and `REAL(h)>0`,
`S_I,pore=-rho*(pi/2)^2*D/max(1e-10,h)*P*f_inhibit`.
Zero P selects this branch, preserving its positive-side pressure derivative;
negative P gets zero diffusion source. This evolves M*P, not P alone.

Dynamic permeability uses
`k=porosity^3*(d_Sauter*sphericity_mean)^2/(150*(1-porosity)^2)`.
For multiple classes, this implicit path evaluates d_Sauter from
`REAL(alpha_s)` and weights sphericity by `REAL(alpha_s)*diam_s^2`.
Sutherland viscosity, when enabled with gas, is
`mu_g=muRef_Suth*(T/Tref_Suth)^1.5*(Tref_Suth+S_mu)/(T+S_mu)` and
`nu_carrier=mu_g/rho_gas`; otherwise the configured carrier viscosity is used.

| Inhibition mode | As implemented in the implicit COMPLEX evaluator |
| --- | --- |
| `OFF` | 1 |
| `STEP` | 1 for `A<maximum_solid_packing`, otherwise 0; real-part comparison |
| `STATIC` | `1-sum_n pascal_coeff(n)*t^(N_inh+n+1)`, `t=clamp((REAL(A)-alpha_trans)/width,0,1)` |
| `DYNAMIC` | Same polynomial, with `alpha_trans_dynamic=10^(-0.28)*REAL(h)^0.12`, capped at `0.99*maximum_solid_packing` if necessary |

`width=maximum_solid_packing-alpha_trans` (or its dynamic counterpart).
The real temporary t and real dynamic threshold discard derivative components
even inside the smooth transition, not merely at its endpoints.
Pascal coefficients are precomputed in the constitutive module during setup.

### Pore coupling outside diffusion

P reduces the basal load only in SI models 1, 11 and 12. It is not subtracted
in the model-9/10 formulas. Explicit bottom/fissural/radial injection updates
M*P using `pore_pres_fract` and the current pressure. Mass exchange transports
P with its mass source, and optional `gas_loss_flag` adds a Darcy-like gas loss
`k/(nu*rho_c)/max(1e-5,h) * (pi/2)*P*f_inhibit`, subsequently limited by
available carrier and packing.

The gas-loss implementation has a solids-presence guard and weights
sphericity by `alpha_s/diam_s`, unlike the implicit diffusion path. Its
DYNAMIC inhibition path also floors h before raising it to 0.12. These are
distinct implementations, not proven interchangeable providers.

## REAL/COMPLEX dependency and branch contract

`eval_jacobian` perturbs only columns selected by `equation_partition%implicit`.
The perturbation is applied to normalized conservative variables:
`q_rel(i)+i*delta`, `q_complex=q_rel_complex*q_org`,
`delta=n_vars*epsilon(wp)`. The matrix is the scaled identity contribution
minus `dt*a_diag*coeff_f*Im(S_I)/delta`. D02 must distinguish derivatives
with respect to q_rel from derivatives with respect to physical q.
The default mask selects both momenta, energy for selector 3, and pore pressure
when enabled; component masses remain explicit. Full physical source
sensitivities to explicit composition/mass are not current Newton columns.

State decoding shares executable include files between REAL and COMPLEX:

| Shared kernel | Active branches / sensitivity contract |
| --- | --- |
| `state_conversion_phys_var.inc` | Wet if `REAL(M)>epsilon`; dry fallback returns h/u/v=0, ambient density/T. For `REAL(M)>eps_sing`, T=E/(M*cp); otherwise ambient T. Nonpositive recovered T also falls back. P=q_pore/M only on the wet path. |
| `velocity_from_conservative.inc` | u=p_x/M above eps_sing; otherwise `sqrt(2)*M*p_x/sqrt(M^4+eps_sing4)` (same for v), retaining the complex perturbation. |
| `mixture_heat_capacity.inc` | Mass-weighted phase heat capacities; carrier fallback below epsilon; gas/liquid flags select the formula. |
| `mixture_from_mass.inc` | Specific-volume mixing, temperature-dependent ideal-gas density and mass-to-volume fractions; carrier fallback below epsilon. |

`complexify` is a complex-step utility, not ordinary complex analysis:
ordering and equality inspect real parts; ABS selects z or -z from REAL(z);
MIN/MAX choose one operand through real-part comparisons and retain the
chosen perturbation. For two-argument MAX, a tie selects the second operand.
Branch-boundary derivatives require a declared one-sided interpretation.
Intrinsic complex SQRT at zero speed and fractional powers near zero require
separate D02 probes, even if a later real-part guard prevents their use.

`c_phys_var` is called inside `IF(rheology_flag)` in the implicit evaluator,
whereas pore diffusion is guarded independently afterwards. The combination
`pore_pressure_flag=T, rheology_flag=F` can therefore read uninitialized local
state there. The pore block also constructs diffusion coefficients before its
wet/pressure guard. These are static control-flow findings requiring targeted
configuration and dry-state tests, not numerical failures demonstrated in D01.

## ODE solutions versus algebraic DIRK solutions

| Term | Current method | What is / is not available |
| --- | --- | --- |
| SI friction | REAL source at the finite-volume candidate; stopping/alignment corrections | Explicitly constructed SI update within the DIRK-weighted stage; not a generic analytic solver |
| Ordinary implicit terms | Newton with complex-step Jacobian | No per-law ANALYTIC dispatch or mixed analytic/Newton groups |
| Model 8 | Closed ODE speed decay over dt | Does not solve `q-r-dt*a_diag*S(q)=0`; slated for removal |
| Model 6 quadratic law | Newton | A collinear flat-bed frozen-property stage could solve `s+k*s^2=r`, but no such production solver exists yet |
| Pore diffusion with frozen coefficients | Newton | An exponential ODE decay differs from the stage relation `P=R/(1+k*dt*a_diag)`; coefficients may depend on state |
| Lava cooling/viscosity | Coupled Newton | T^4 cooling and T-dependent momentum prevent inferring a generic closed update from a scalar special case |

For the model-6 illustrative scalar relation, `s` and `r` are momentum
magnitudes and `k=dt*a_diag*friction_factor/(rho*h^2)` only under ordinary
velocity recovery, flat bed and fixed rho/h. This is a candidate derivation,
not an implemented method or acceptance of D13. General geometry, near-dry
regularization, coupling, positivity and group isolation must be established.

The current stage path is: build finite-volume candidate; evaluate REAL SI;
clip sign reversals/align momentum and reconstruct effective SI increment;
run Newton; apply further stopping/alignment corrections; store the **effective**
implicit increment `(q_guess-q_si)/(dt*a_diag)`. Newton differentiates the
implicit source only, not this whole corrected map. D02 source tests and later
D09–D14 stage tests must remain distinct.

## Findings and D02 worklist

No item below is fixed or numerically certified by this documentation change.

| ID | Static observation / risk | Next evidence or action |
| --- | --- | --- |
| D01-F01 | Model 7 is absent and rejected; no identity is established | Keep unsupported in D08 mapping; add parser negative coverage when changing input handling |
| D01-F02 | Model 8 integrates an ODE using dt without a_diag | D04 removal, including solver/integrator/parser/shear diagnostics; D13 is separate |
| D01-F03 | `xi_temp` is shared mutable scratch in models 1/12; `T_env=300` is a shared write in model 3 | D03 ownership/thread audit; D02 repeated calls and requested versus actual thermal settings |
| D01-F04 | STATIC/DYNAMIC inhibition and multi-class permeability project intermediate quantities to REAL | Compare full physical finite differences against CS; separate current Newton columns from lost explicit-column sensitivities |
| D01-F05 | Diffusion lacks the gas-loss solids guard; sphericity weights differ between the two paths | Multi-class, zero-solid and packing-limit fixtures; do not merge formulas without a scientific decision |
| D01-F06 | Pore-only configuration can skip physical-state initialization; dry diffusion computes coefficients before its guard | Targeted pore-only/dry calls with strict arithmetic/runtime checks; define supported flag combinations |
| D01-F07 | Several validations use AND for missing/negative parameters (1, 9, 10); xi/Fr_0/U_w/I_0 can reach problematic zero limits | Negative/partial/zero input matrix; document expected limits before changing validation |
| D01-F08 | Model-9 comments reverse the formula's endpoint values; models 9/10 do not compute curvature locally and differ in g' guards | Independent endpoint/load tests; retain current behavior until a separately validated correction |
| D01-F09 | Thresholds and floors differ across implicit/SI laws; zero-speed SQRT and h powers are piecewise | One-sided threshold probes, wet/dry limits, pure x/y and rotated velocity states |
| D01-F10 | Thermal input surface exceeds active cooling terms; lahar viscosity is duplicated; diagnostics have their own load formulas | Trace D06 clients; source/diagnostic distinction; avoid claiming dormant conduction/heating support |

D02 minimum matrix, to be implemented and accepted separately:

1. Every supported law: regular moving state, rest, h approaching zero, both
   axes and oblique flow, slope on/off; curvature and P where actually used.
2. Implicit models 1–6/12 and pore diffusion: REAL versus real part of COMPLEX,
   zero imaginary residual at unperturbed inputs, source derivatives and scaled
   Newton columns, finite-difference step sweeps away from branches. Do not
   report zero implicit model-9 output as a derivative test of mu(Fr).
3. SI-only terms: independent REAL formulas and limiting/dissipation checks;
   no imaginary-perturbation assertion through a nonexistent COMPLEX API.
4. Lava: temperature/reference/ambient values, radiation and convection
   separately, h=1e-10 neighborhood. Lahar: Tc=20 from both sides, solid
   concentration, thresholds and deposition-viscosity consistency.
5. Pore: P negative/zero/positive, h=1e-10 neighborhood, all four inhibition
   modes, static/dynamic permeability, mono/multi-class solids, A approaching
   0 and packing, Sutherland on/off, gas-loss interactions and flag validity.
6. Stochastic inputs: frozen Z for differentiation, xi clamp endpoints,
   Fr shift/clamp and h=0.01, Lucas U_w and h=0.1; RNG evolution is a separate
   stochastic-module contract. No random draws belong inside source probes.
7. Integration checks for source stopping/alignment and the stored effective
   RK increment after corrections, distinguished from raw source Jacobians.

Existing reusable coverage is limited: `TEST_RHEOLOGY_FROUDE` checks model 9
domain/REAL formula and zero implicit REAL/COMPLEX output;
`TEST_STATE_CONVERSION` checks shared conversion paths and some derivatives;
`TEST_IMEX_STAGES` includes a model-6 implicit-drag fixture;
lava setup/example and source-mass tests exercise integration/configuration.
None constitutes complete per-law D02 characterization. Existing tests were
inspected, **not rerun** as part of D01.

## Validation and acceptance boundary

Classification: documentation/static audit, T1. Checks validate source-hash
provenance, routine/range/dispatch references, parameter declarations, local
document links, JSON structure, model coverage and whitespace. Numerical T2–T4,
builds, CI and simulations are excluded because no executable or input changes.
The report's assertions concern the pinned revision, not future HEADs.

To recheck the pinned file contents from the repository root (a later source
change is expected to invalidate the corresponding hashes):

```sh
python3 -B - <<'PY'
import hashlib
import json
from pathlib import Path

path = Path('docs/acceptance/results_d01_inventory_2026_10_21.json')
inventory = json.loads(path.read_text())
for item in inventory['source_files']:
    actual = hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()
    assert actual == item['sha256'], item['path']
print('PASS: D01 source hashes')
PY
git diff --check
```

D01 is complete for the formula/branch/dependency inventory. D02 remains
**DA FARE**; numerical REAL/COMPLEX parity, derivatives, physical acceptance
and later extraction equivalence are not established by this static gate.
D03 remains the full symbol-level ownership audit. No commit is made by this task.
