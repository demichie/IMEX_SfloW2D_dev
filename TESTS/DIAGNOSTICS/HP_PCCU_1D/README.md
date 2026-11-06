# HP-PCCU 1-D numerical diagnostics

These files support reproducible investigations of hydrostatic-path PCCU
spatial discretization, shallow-flow reconstruction and wet/dry behavior
on simplified one-dimensional tests.

| File | Purpose |
| --- | --- |
| `hp_pccu_1d_core.py` | Python reference implementation of the baseline hydrostatic-path PCCU operator. |
| `kt_v2_1d.py` | Reduced pressure-split 1-D diagnostic with independently reconstructed transport variables. |
| `run_inclined_1d.py` | Inclined-excavation experiment driven by the reference Python core. |
| `smooth_bump_reconstruction_1d.py` | Smooth-bump reference/reconstruction comparison. |
| [`fortran_inclined_1d/`](fortran_inclined_1d/README.md) | Matching input/DEM/initial-state fixture for the compiled Fortran solver. |

All Python commands below are run from this directory and require
Python 3 and NumPy. Optional plotting may need Matplotlib.

Run the inclined-excavation Python diagnostic:

```bash
python3 run_inclined_1d.py --help
```

Run the reduced KT/HP-PCCU diagnostic:

```bash
python3 kt_v2_1d.py --help
```

Run the smooth-bump reconstruction diagnostic:

```bash
python3 smooth_bump_reconstruction_1d.py --help
```

The nested Fortran fixture documents its own preparation and execution
commands. The algorithms here were **relocated, not refactored**.
