# D17.9 — Explicit stage routing

Status: **implemented; local GNU validation passed**. This implements roadmap
D17.9 and D17.19 implementation item 11. It routes the existing explicit cell
rates and boundary-state hook through the validated `stage`/`explicit` group;
it does not change the physics, RK weights, time argument, source windows or
spatial flux formulation.

## Implementation

- `eval_expl_terms` now requires an explicit stage group and gates curvature,
  circular bottom injection, fissural injection and radial cell injection by
  their registered `term_id`s. The existing physical activation flags, cell
  coverage, wet/dry checks, pulse windows and rate assembly remain in force.
- `time_integration_2d` selects all explicit stage groups once before the cell
  loops, evaluates their cell rates in registry order, and sums the returned
  rates into the current EX stage. No group is inferred from equation-array
  position.
- Boundary injection remains a primitive boundary state consumed by
  reconstruction and the spatial flux, not a cell rate. Its membership in the
  same explicit stage group now gates both CFL reconstruction and RK spatial
  evaluation. Legacy direct reconstruction calls that omit the optional gate
  retain their prior flag-based behavior.
- `step` mass-exchange groups and their common limiter/application batch remain
  on the existing post-RK path; they are not evaluated by the stage interface.

## Local validation

| Check | Result |
|---|---|
| Curvature/source regression | Strict Fortran test passes; selected curvature remains active and an unselected circular injection is excluded |
| D17 source registry | Strict and optimized: 864 legacy configurations, 14,857 assertions, 28 invalid configurations rejected per profile |
| Spatial operator | Strict and optimized; all existing topography/CFL reference checks pass at 1 and 4 threads |
| IMEX stages | Strict and optimized; raw-stage and CFL regressions, exact-flow and mixed ANALYTIC/NEWTON trajectories for `N_RK=2/3/4` pass |
| Boundary reconstruction | Strict and optimized; selected and unselected source traces pass for W/E/S/N and collapsed-grid cases at 1 and 4 threads |
| Production executable | Out-of-tree GNU Fortran build and link pass with the installed NetCDF libraries |
| Whitespace | `git -c core.whitespace=cr-at-eol diff --check` passes |

No Intel ifx, Linux, CI or cross-platform acceptance is claimed. Subsequent
group-filtered semi-implicit/implicit routing is documented in
[D17.10 acceptance](D17_10_GROUP_FILTERED_IMPLICIT_RATES.md). Later roadmap
requirements remain open; this report does not mark D17 complete.
