# Etna lava example preparation and validation

`EXAMPLE_ETNA_LAVA` replaces `EXAMPLE_LAVA_IMO` as the lava case selected by
the user for N9; the other selected long case remains `EXAMPLE_ETNA`.
The IMO example and its historical evidence are retained. N7/N8 remain closed;
this package does not approve the N9 scientific baseline.

## Input and original topography

The [example](../../EXAMPLES/EXAMPLE_ETNA_LAVA/README.md) preserves the supplied
input byte for byte, including its original output flags, 3600 s duration,
60 s output interval, 250 by 150 cells at 10 m, rheology 3 and two fissures.
Slope correction and curvature remain OFF. Omitted `CFL_POLICY` selects
`directional_min`; the unchanged input uses `CFL=0.24` and `N_RK=2`.

`prepare.py` streams the original raster from the already tracked
`EXAMPLE_ETNA/DEM/Etna2014_crop.zip`. The extracted bytes match the original
raster in the user-provided lava example. No ellipsoid/niche preprocessing,
smoothing or initial-volume generation is performed. The solver retains its
normal topography import. The simulation extent is inside the original 3000
by 2000, 1 m DEM. [fixture.json](../../EXAMPLES/EXAMPLE_ETNA_LAVA/fixture.json)
pins input, archive and raster SHA-256 values.

Preparation is idempotent and refuses to overwrite differing input/DEM data.
The archive is shared rather than duplicated; the generated 102 MB raster and
simulation outputs are ignored locally by Git. No existing user outputs are
copied, overwritten or removed.

## Full duration execution checks

The [execution evidence](results_etna_lava_example_2026_10_11.json) records four
complete runs to 3600 s: strict/debug and optimized builds, each with actual
one/four-thread teams. All four pass the example validation checks. The
executables are those from the completed CFL audit; their production-source
hashes were rechecked against commit `8fc540eb` before execution. This package
does not claim another fresh compilation.

Every run preserves the exact input and original DEM hashes and produces:

- All 61 NetCDF times from 0 to 3600 s and 61 thickness/61 temperature ESRI rasters.
- Finite decoded numeric variables without NetCDF fill values, nonnegative
  thickness, expected coordinates/shapes/units, unchanged bed and consistent
  `eta=b+h` within the existing 8192-epsilon scaled roundoff criterion.
- A final binary checkpoint. Within each compiler profile, all decoded NetCDF
  variables and checkpoint bytes are exactly identical across one/four threads.

The eight preparation checks cover real-raster extraction, input identity,
domain coverage, idempotence, bad checksums, path traversal, symbolic links and
preservation of differing existing files. All 16 acceptance-tool checks pass.
No numerical source, reference kernel or existing acceptance tolerance changes.

## Source budget and remaining N9 work

The two rectangular fissures have nominal areas 600 and 500 m2. At the supplied
0.0084 m/s injection velocity, their nominal total supply is 9.24 m3/s, or
33264 m3 over the full run. The final volume retained in the domain is about
32734.01306 m3, a difference of 529.98694 m3. Some boundary cells are wet
(maximum recorded boundary thickness about 0.17085 m).

This is not a verified closed-domain mass defect: the N9 budget must include
the actual integrated numerical boundary flux. That flux has not yet been
integrated, so outflow remains a hypothesis requiring the balance check, not
a declared explanation or an approved budget tolerance. No domain extension,
shortened duration, parameter change or tolerance relaxation is used here.

The subsequent [N9 long-case checks](N9_LONG_CASES_2026_10_11.md) include ETNA and
actual full-duration continuation. The user limits the lava budget to the last
saved output with dry boundary cells (3180 s), before the first positive-boundary
output at 3240 s. That precontact volume budget agrees at roundoff scale; the
postcontact final-domain difference is outside the authorized mass-check scope.
Physical fields, ESRI rasters and checkpoints match uninterrupted runs exactly,
but diagnostic CSV history fails on continuation. N9 remains open for the
diagnostic correction and explicit scientific baseline approval.
