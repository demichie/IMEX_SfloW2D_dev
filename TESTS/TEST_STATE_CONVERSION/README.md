# State-conversion characterization test

This test verifies the state-conversion routines after the constitutive split,
the introduction of the shared REAL/COMPLEX core, and the correction of the
mixture heat-capacity weighting.

The passing baseline covers:

- `q -> qp -> q` round trips with `alpha` and `h*alpha` storage;
- energy and temperature-transport configurations;
- gas, two solid classes, one additional gas, liquid, stochastic transport,
  and pore pressure;
- dry and near-dry states;
- exact gas-liquid-solid energy and thermodynamic round trips;
- REAL/COMPLEX values at an unperturbed wet state;
- complex-step derivatives compared with centered finite differences.

The heat-capacity correction removes the previous gas-liquid-solid discrepancy:
each carrier-component mass fraction now contributes exactly once to
`c_p,mix`, and the corresponding `q -> qp -> q` and REAL/COMPLEX value checks
must agree to round-off accuracy.

The shared conversion formulas use explicit element-wise multiplication and
`SUM` for composition-weighted quantities. Unlike Fortran `DOT_PRODUCT`, this
does not conjugate a complex first argument, so the COMPLEX implementation is
an analytic continuation of the REAL formulas. All tested complex-step
derivatives must now agree with centered finite differences within the stated
tolerance.
