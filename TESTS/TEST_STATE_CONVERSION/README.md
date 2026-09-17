# State-conversion characterization test

This test freezes the behavior of the state-conversion routines before the
`constitutive_2d` module is split. It does not correct constitutive formulas.

The passing baseline covers:

- `q -> qp -> q` round trips with `alpha` and `h*alpha` storage;
- energy and temperature-transport configurations;
- gas, two solid classes, one additional gas, liquid, stochastic transport,
  and pore pressure;
- dry and near-dry states;
- REAL/COMPLEX values at an unperturbed wet state;
- complex-step derivatives compared with centered finite differences.

Two existing differences are intentional expectations of this baseline:

1. Gas-liquid-solid conversion does not exactly round-trip the energy and
   thermodynamic fields.
2. Nine tested REAL/COMPLEX derivatives differ, although the corresponding
   unperturbed gas-solid values agree.

The test reports these as `KNOWN DISCREPANCY` and requires their count and
magnitude to remain bounded. M12 can then change these expectations in the
same commit that introduces the shared REAL/COMPLEX conversion core.
