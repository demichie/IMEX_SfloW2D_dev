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

One independent known difference remains: nine tested REAL/COMPLEX derivatives
differ, although the corresponding unperturbed values agree. The cause is the
implicit conjugation performed by Fortran `DOT_PRODUCT` when its first argument
is complex; correcting those products belongs to a separate test-backed commit.

The test reports those derivative differences as `KNOWN DISCREPANCY` and
requires their count and magnitude to remain bounded until that correction is
applied.
