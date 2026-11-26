# Two-layer hyperbolicity: characteristic limits and acceptance

**Status (26 November 2026):** the present formulation is accepted as an
**explicitly limited experimental/research model**, not as a generally
hyperbolic or fully physically validated model. This document makes no
change to the governing equations, CFL, sources or validation gate.

## 1. Scope and implemented operators

The following derivation is based on the actual Fortran implementation:

- [equation_terms_2d.f90](../src/physics/equation_terms_2d.f90):
  inertial flux `F_n=u_n*q`, density and reduced gravity.
- [pccu_2d.f90](../src/spatial/pccu_2d.f90):
  hydrostatic one-form `-G*(Gamma*h*d(eta)+0.5*h**2*d(Gamma))`;
  interlayer products `-G*h1*d(L)` and `-L*d(G*h1)`.
- [hyperbolic_2d.f90](../src/spatial/hyperbolic_2d.f90):
  interlayer pressure corrections and total-depth face-speed bounds.
- [spatial_operator_2d.f90](../src/spatial/spatial_operator_2d.f90):
  spatial residual assembly for the two layers.

Here `L=Gamma2*h2`, `Gamma_i=rho_i*gprime_i` and
`gprime_i=g*(1-rho_a/rho_i)`, where `rho_a` is atmospheric density.
The reconstructed surfaces are `eta1=B+h1`, `eta2=B+h1+h2`.
The geometry factor is `G=1/(1+B_x**2+B_y**2)`.

Analyze a **smooth, locally planar, homogeneous wet/wet state**:
`h1,h2>0`, `rho1,rho2>rho_a`, a constant positive `G`,
and locally frozen temperature/composition in each layer.
The L2 nonnegative-pressure clamp is therefore inactive.
Choose a normal direction with layer velocities `U1,U2`.
This is a local continuous principal-symbol calculation, not an
eigenanalysis of the discrete PCCU flux or a certification of the
full thermodynamic/wet-dry system.

The homogeneous depth/normal-velocity subsystem is

```text
d_t h1 + U1*d_x h1 + h1*d_x U1 = 0
d_t h2 + U2*d_x h2 + h2*d_x U2 = 0
d_t U1 + U1*d_x U1 + G*gprime1*d_x h1
                         + G*(Gamma2/rho1)*d_x h2 = 0
d_t U2 + U2*d_x U2 + G*gprime2*d_x(h1+h2) = 0
```

With `W=(h1,h2,U1,U2)`, the gravitational principal block is

```text
       [ U1          0          h1  0  ]
A_n =  [ 0           U2         0   h2 ]
       [ G*gprime1   G*Gamma2/rho1 U1 0 ]
       [ G*gprime2   G*gprime2   0   U2 ]
```

Composition, thermal and tangential-velocity transport adds advected
modes and composition-gradient pressure couplings; sources affect
the evolving local state, but not this frozen-state gravitational block.

## 2. Gravitational eigenvalues and density ordering

Set `a=G*gprime1*h1`, `d=G*gprime2*h2`,
`k=G**2*(Gamma2/rho1)*gprime2*h1*h2`.
The **exact quartic of this gravitational block** is

```text
P(lambda) = ((lambda-U1)**2-a)*((lambda-U2)**2-d)-k = 0.
```

When `U1=U2=U`:

```text
mu_+/- = (a+d +/- sqrt((a-d)**2+4*k))/2
lambda = U +/- sqrt(mu_+), U +/- sqrt(mu_-)

a*d-k = G**2*gprime2*h1*h2*g*(rho1-rho2)/rho1.
```

Under the above positive-depth/positive-reduced-gravity assumptions:

- **`rho1>rho2>rho_a`:** the four gravitational eigenvalues are
  real and distinct at zero slip (stable density ordering).
- **`rho1=rho2>rho_a`:** `mu_-=0`; `lambda=U` is double.
  For two positive depths the gravitational block is defective,
  hence it is **not strongly hyperbolic**, even though its
  eigenvalues are real.
- **`rho2>rho1>rho_a`:** `mu_-<0`; one complex-conjugate
  pair exists already at zero slip (inverted density ordering).

These statements concern the gravitational subsystem, not a
claim about cells visited by any particular simulation.

## 3. Interlayer shear and complex roots

With unequal normal velocities, even the stably stratified case
may have complex eigenvalues. Let `DeltaU=U2-U1`. At a transition
where two real roots merge, the **exact implicit threshold** satisfies

```text
P(lambda_*, DeltaU_crit) = 0
dP/dlambda(lambda_*, DeltaU_crit) = 0.
```

The first positive `abs(DeltaU)` where all four roots cease to be
real and distinct depends on depths, densities, gravity and the
geometric factor. It is not a model-wide constant.

For a *small* relative density contrast the leading-order estimate is

```text
DeltaU_crit**2 ~= G*g*(rho1-rho2)/rho1*(h1+h2).
```

Do **not** replace the exact quartic by this asymptotic formula outside
its regime of validity.

Example: `G=1`, `g=9.81 m/s**2`, `rho_a=1.2 kg/m**3`,
`h1=h2=1 m`, `rho1=1000 kg/m**3`, `rho2=200 kg/m**3`.
Numerically solving `P=P'=0` gives two positive slip transitions,
approximately **4.6535** and **7.5194 m/s**. In between, a
complex-conjugate pair occurs in the chosen 1D direction.
Reappearance of real roots for larger 1D normal slip does **not**
certify two-dimensional hyperbolicity.

For a flat/isotropic local state, the normal slip for a horizontal
unit direction `n=(nx,ny)` is

```text
DeltaU_n = (u2-u1)*nx + (v2-v1)*ny.
```

A 2D hyperbolicity claim requires considering **all** directions.
Even if a large vector slip gives real roots along its own direction,
another oblique direction can project it into the complex-root band.
The different slope metrics and dry/wet faces are not covered by
this homogeneous flat-state directional argument.

## 4. The total-depth CFL envelope is conditional

The solver calls
[eval_total_depth_face_speed_bounds](../src/spatial/hyperbolic_2d.f90)
and forms signed bounds using the present-layer speeds evaluated
at `H=h1+h2`:

```text
Ui +/- sqrt(G*gprime_i*H)
a_minus = min(0, all face-endpoint/layer lower bounds)
a_plus  = max(0, all face-endpoint/layer upper bounds).
```

For the planar wet/wet block with `rho1>=rho2>rho_a`,
**every real gravitational eigenvalue** satisfies

```text
min_i(Ui-sqrt(G*gprime_i*H))
    <= lambda <=
max_i(Ui+sqrt(G*gprime_i*H)).
```

Indeed, above the upper envelope, each factor
`(lambda-Ui)**2-G*gprime_i*hi` exceeds
`G*gprime_i*h_other`, and their product exceeds `k` since
`Gamma2/rho1<=gprime1` for `rho2<=rho1`.
Thus `P(lambda)>0`; the lower case is symmetric.
This justifies the **containment of real roots in this restricted
regime**, not the reality of all roots, strong hyperbolicity, or
stability of the full discrete method. No corresponding general
claim is made for inverted ordering or nonpositive reduced gravity.

**Reducing the timestep cannot make complex eigenvalues of the
continuous principal symbol real.** Passing a CFL check does not
certify hyperbolicity.

## 5. Model interpretation and explicit limitations

- Both depths must be positive for this derivation; drying, emergence
  of L2, and nonpositive reduced gravity require separate analysis.
- The L2 buoyancy clamp, variable composition/temperature gradients,
  full convective multiplicities, spatial reconstruction and
  nonuniform terrain are outside this specific proof.
- Elutriation using the dynamic L1 donor composition may form an
  L2 with density close to L1. The transfer also carries donor
  momentum. Atmospheric entrainment can subsequently dilute L2,
  but neither sufficiently fast dilution nor small ongoing shear
  is guaranteed by the equations.
- Moving interlayer exchange before deposition/entrainment changes
  the **source trajectory**, not the continuous characteristic
  polynomial. No Lie ordering constitutes a hyperbolicity repair.
- Existing positivity, conservation, restart, OpenMP and integration
  tests do not certify real characteristic roots in every cell.
  No verified full-field spectral classification of the Etna
  simulation is available at this decision point.

## 6. Limited scientific acceptance

**Accepted for experimental use:** the current two-layer equations,
the existing total-depth CFL estimate, both implemented elutriation
composition options (`-1`: local L1 donor fraction; positive:
prescribed lifted fraction), and the existing source-group
configuration and physical closures, with the explicit restrictions
and uncertainties above.

**Not established:** general strong hyperbolicity for arbitrary
admissible states, validity of inverted-density or excessive-shear
regimes, physical reliability of solutions in complex-eigenvalue
regions, and quantitative validation of all interlayer closures.

**No change to the physical or numerical model is authorized by this
decision:** do not add ad hoc shear/density thresholds, artificial
pressure corrections, numerical hyperbolicity repairs or new scientific
milestones. Both `N_LAYERS=1` and `N_LAYERS=2` are directly selectable
in the input file; no special environment activation is needed.
TL13a is functionally complete; direct execution of two-layer cases
does **not** imply general TL13b scientific approval.

For interpretation of an individual run, the wet/wet density ordering
and its characteristic roots should be assessed using its actual
layer-resolved states rather than inferred from successful completion.

See [TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
and [TL13B_C_PHYSICAL_CLOSURES.md](TL13B_C_PHYSICAL_CLOSURES.md).
