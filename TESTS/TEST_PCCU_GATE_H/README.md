# HP-PCCU Gate-H refinement regression

This test reproduces the compact four-sided one-cell inclined excavation from
the frozen HP-PCCU validation package on the three square grids
`dx = dy = 0.4, 0.2, 0.1 m`. The nominal 8 m by 6 m excavation is translated
by half a cell in both directions and constructed from one globally continuous
Q1 vertex field. Its initial thickness is the exact Q1 cell average, so the
initial free surface is the unexcavated inclined plane.

Two layers are exercised:

- baseline HP-PCCU with `G=1`;
- slope correction and the separate curvature source both enabled.

`reference_metrics.csv` contains field fingerprints generated with
`HP_PCCU_clean_reference_cores_2026-09-26_drySafe_rebase.zip`. The test checks
positivity, conservation, zero uphill transfer, the G=1 lateral symmetry, mass
partitions, thickness moments, total variation and extrema. Aggregate
fingerprints keep the regression small while making it sensitive to changes in
the complete 2-D thickness field.

The reference core SHA-256 hashes, as recorded in the package, are:

- `ffcb695ff8d8d2bd06e433d70b2a519c0b0012d959eb664b617deea83e8f64e2`
  for `hp_pccu_1d_core.py`;
- `c63331366cadf9f04d802a3f058088b9c175399b2d6dc7015068a0ad110f1577`
  for `hp_pccu_2d_core.py`.

During the port audit the complete Fortran and Python arrays were also compared.
For slope plus curvature, the normalized errors were:

| dx | L1(h) | L1(eta) | Linf(h) [m] |
|---:|---:|---:|---:|
| 0.4 | 7.1697e-5 | 2.7927e-5 | 1.0116e-3 |
| 0.2 | 6.5239e-5 | 2.5412e-5 | 8.6013e-4 |
| 0.1 | 6.4436e-5 | 2.5099e-5 | 1.2830e-3 |

Run from the repository root with:

```sh
TESTS/TEST_PCCU_GATE_H/run_test.sh /path/to/IMEX_SfloW2D
```
