# EXAMPLES README migration — second-pass assessment (2026-10-16)

Source repository: `demichie/IMEX_SfloW2D_dev` at `07b740f7`.

The Examples tree contains 20 standardized `README.md` files and 19 legacy `README.txt` files. `EXAMPLE_ETNA_LAVA` has no legacy README; no Markdown link points to a missing `README.txt`. The legacy files, scripts, input namelists and DEM data remain preserved.

## Validation carried out

- Checked the inventory, YAML metadata, required headings, balanced Markdown code fences and legacy README links for all 20 Markdown files.
- Regenerated the public example catalogue from this tree: 20 pages, zero content warnings; MkDocs `--strict` completed successfully with the generated configuration.
- Reviewed parameter semantics against generators, inputs and solver code for the complex source cases below. This is not a lossless, section-by-section comparison of all 19 legacy texts.
- From a fresh extraction at `07b740f7`, `EXAMPLE_2D`, `EXAMPLE_BW2D`, `EXAMPLE_ENGWELL2016`, `EXAMPLE_TAAL` and `EXAMPLE_TAAL_BOTTOM` reached their configured end times and wrote nonempty NetCDF output. The seven previously repaired NumPy 2.x generators and the Etna ensemble were validated in earlier C06 runs.
- Runtime exceptions are recorded below: the documented `EXAMPLE_DAMBREAK` pore-pressure case fails, the documented `EXAMPLE_LATERAL_SOURCE` input produces a non-finite source and zero mass output despite exit status 0, and `EXAMPLE_CF` did not complete within the run window.

## Corrections and caveats

- `EXAMPLE_BW2D`: eight script arguments; conflicting pulse narratives not silently erased.
- `EXAMPLE_LATERAL_SOURCE`: the third argument is `MFR_SOURCE` in kg/s and the sixth is written to `XS_SOURCE` as a mass fraction. The documented combination (`T_SOURCE=900 K`, `XS_SOURCE=0.01`) gives source reduced gravity `-19.33 m/s²`, then NaN source thickness/velocity/mass flow and zero simulated mass; the solver nevertheless exits 0. A diagnostic run with `XS_SOURCE=0.8` made the derived source finite but still produced zero mass, so it is not a validated replacement configuration.
- `EXAMPLE_ENGWELL2016`: the sixth generator argument is written to `XS_SOURCE` as a mass fraction; the legacy volume-fraction description is incorrect.
- `EXAMPLE_DAMBREAK`: current generator is `create_dambreak.py`; historical filenames are identified as such. The documented `PORE_PRESSURE_FRACT=1.0` run reaches NaNs and fails with `negative thickness after IMEX Runge-Kutta update` at about 0.36 s. A 0.5 fraction fails similarly at about 0.49 s; a 0.0 fraction reached the configured 3 s, but removes the initial excess pore pressure and is not equivalent to the documented experiment. With fraction 1.0, the generated pore pressure exactly cancels the initial `rho_m * g_red * h` vertical stress; the solver clamps effective stress to zero. Disabling dynamic permeability or gas loss separately, and a diagnostic switch to rheology model 6 with a valid friction factor, still fail at about 0.36 s. This points to the pressure/effective-stress path rather than either permeability option or rheology model 11 alone; the precise numerical failure mechanism remains unresolved. The generator prints the full hydrostatic excess pressure before multiplying by the requested fraction, so its diagnostic print does not show the initialized value when the fraction is below 1.
- `EXAMPLE_ETNA`: `cleanFolder.sh` capitalization and ellipsoid-preprocessing explanation.
- `EXAMPLE_2D`: old `create_example4.py` is not the committed generator filename. The current generator needs NumPy and only imports Matplotlib when plotting is requested; it does not use the Basemap dependency listed by the legacy README. Its generated input enables both NetCDF and ESRI output.
- `EXAMPLE_ENSEMBLE_ETNA`: placeholder `ENSEMBLE_`, link/template execution stages retained.
- `EXAMPLE_CF`, `EXAMPLE_TAAL` and `EXAMPLE_TAAL_BOTTOM`: their committed namelists have ESRI output disabled, contrary to the generic legacy output descriptions.
- `EXAMPLE_CF`: the extracted run reached a displayed time of 105.584 s, with `dt` rounded to 0.00000 and no visible time advance during the final several minutes. It was stopped after about 9 minutes; the run did not reach its configured end time.
- `EXAMPLE_TAAL_BOTTOM`: because `BOTTOM_RADIAL_SOURCE_FLAG=T`, `TIME_PARAM(4)` is a phase offset, not a total source duration or cut-off. Its current value (150 s) is five 30 s periods and therefore gives no net phase shift. This differs from lateral `EXAMPLE_TAAL`, where `TIME_PARAM(4)` is an absolute cut-off.
- `EXAMPLE_ETNA_LAVA`: extra DEM and acceptance context restored; validation relates to a pinned source and input configuration.

## Outstanding work

This pass **does not certify lossless transfer of every original sentence**. Keep the 19 legacy files until a section-level completeness comparison of all paragraphs, parameter tables, references and usage commands is completed. The default lateral-source recipe and the dam-break pore-pressure case need a scientifically approved correction or explicit acceptance as known limitations; the CF runtime needs a successful completion or a documented, reproducible explanation for its timestep collapse. The remaining examples and optional postprocessing dependencies also need validation. Do not claim that all 20 solver runs have passed.

The MkDocs generator should render the entire Markdown body rather than extracting only required sections. No legacy README files were removed.
