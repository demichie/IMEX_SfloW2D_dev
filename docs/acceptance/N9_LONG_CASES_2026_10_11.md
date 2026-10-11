# N9 long cases and restart diagnostics

> Historical v4 characterization: this report preserves the diagnostic
> continuation failures observed before checkpoint v5. The defects were
> subsequently corrected and revalidated; see
> [the v5 restart report](N9_RESTART_V5_2026_10_11.md) and
> [v5 evidence](results_n9_restart_v5_2026_10_11.json). Do not interpret the
> failures below as the current v5 result.

## Historical v4 results (superseded for restart diagnostics)

The selected ETNA and ETNA_LAVA cases preserve their full durations, physical
inputs and default `directional_min` CFL policy. All eight physical restart
comparisons pass, but the runout and mass-centre CSV continuation checks fail.
**N9 remains open**: those diagnostic defects must be resolved before scientific
baseline approval and semantic refactoring.

## Source and comparison scope

The source revision is `f6a004fa`, with unchanged production sources from the
audited `8fc540eb` numerical implementation. The two executables are reused
from its clean strict/debug and historical optimized builds, not relabeled as
fresh builds. Every production source hash is checked against both the build
manifest and the preceding lava evidence before execution or characterization.
The existing 45/45 CFL acceptance audit and N7/N8 closure remain applicable;
that full audit was not rerun for this documentation and observation package.

The original ETNA preprocessing script and parameters produce `new_dem.asc`,
`init_volume.asc` and `ellipsoid.asc` byte-identical to the user's existing
rasters. A test-only NumPy ASCII table reader replaces the script's sole pandas
read, leaving geometric expressions and output formatting unchanged. Input,
scripts, archive, raw raster and generated-raster hashes are recorded.
ETNA_LAVA retains the supplied input and unmodified original DEM.

The September 26 dry-safe reference archive matches both clean-core hashes in
the frozen Gate-H contract and repository. The earlier September 24 archive is
not substituted for that reference. This provenance check does not rerun the
external package's complete scientific suite.

## Mass budgets

ETNA has no wet boundary cells in the saved outputs through 200 s. Its initial
positive depths exceed the dry cutoff. Initial native fields therefore contain
the entire initial mass, and the final budget uses all conservative checkpoint
cells, including positive dry residues. Total and solid-mass endpoint differences
are at accumulated floating-point rounding scale.

| Case and build | Budget interval [s] | Total mass residual [kg] |
| --- | --- | ---: |
| ETNA optimized | 0 to 200, initial versus conservative endpoint | -1.4901161e-8 |
| ETNA strict/debug | 0 to 200, initial versus conservative endpoint | -1.7434359e-6 |
| ETNA_LAVA optimized | Saved outputs 0 to 3180, maximum absolute residual | 8.9406967e-8 |
| ETNA_LAVA strict/debug | Saved outputs 0 to 3180, maximum absolute residual | 2.8312206e-7 |

For ETNA, native `rhomix` and solid fractions intentionally vanish below the
existing dry-depth cutoff, while conservative mass is retained. Consequently,
intermediate sums of `h*rhomix` omit inventories of fractions of a gram and are
recorded as masked-output observations, not conservative budget failures or
proof of full conservative conservation at every intermediate stage.

For the pure liquid ETNA_LAVA case, density is constant at 2600 kg/m3. Using
`2600*sum(h)*cell_area` retains dry residuals too. The two fissures supply a
nominal 9.24 m3/s; at 3180 s the expected volume is 29383.2 m3. Maximum absolute
precontact volume residuals are 3.2741809e-11 m3 (optimized) and 1.0913936e-10 m3
(strict/debug). Results are identical for one and four threads within each build.

The user authorizes the lava budget only before the flow reaches boundary cells.
3180 s is the last saved output with exactly zero boundary depths; 3240 s is the
first saved output with positive boundary depths. The output cadence is 60 s,
so the exact first-contact time is not determined. The approximately 530 m3
final-domain volume difference at 3600 s is outside this budget scope; it is not
declared a mass defect or a proven outflow balance. Full-duration field and
restart checks still extend to 3600 s. No scientific tolerance is relaxed or
new mass threshold inferred from these observations.

## Physical continuation and output schema

Both builds use actual OpenMP teams of one and four threads. ETNA is split at
the existing 100 s output event and continued to 200 s; ETNA_LAVA is split at
1800 s and continued to 3600 s. Temporary inputs change only first-leg `T_END`
and continuation `RESTART`/`RESTART_FILES`; physical and numerical parameters,
output flags and output intervals are unchanged.

All eight comparisons have bitwise-identical decoded NetCDF arrays, identical
metadata, identical ESRI rasters and byte-identical final version-4 checkpoints.
All four full/resumed one-versus-four-thread pairs also pass for fields and
checkpoints. Native outputs have 41 ETNA snapshots and 61 ETNA_LAVA snapshots,
with 24 and 20 numeric variables respectively. Fields are finite, thickness is
nonnegative, the bed is unchanged, and coordinate/time arrays match the inputs.

The native schema distinguishes free-surface elevation `eta` from vertical
velocity `w`; units, names, rank and dimension ordering are checked. The decoded
checkpoint is scoped to GNU Fortran little-endian sequential records, validates
magic/version/layout/model flags and checks conservative mass, component mass
and thermal-energy identities against the physical output. Dry diagnostic
density/fractions are checked under their existing zero-output convention rather
than incorrectly demanding a conservative identity from masked fields.
The field-identity checks retain the existing 8192-epsilon roundoff scale;
continuation and thread comparisons allow no field differences. Cross-build
bitwise equivalence and universal compiler portability are not claimed.

## Failed diagnostic continuation

All eight CSV continuation checks fail without changing limits:

- ETNA loses its initial 0 to 100 s CSV history and the headers; resumed files
  contain only 101 to 200 s. Runout differs by up to about 4206686 m at common
  times, and directional mass-centre extents use the wrong reference.
- ETNA_LAVA loses its initial 0 to 1800 s CSV history and the headers. Values
  at the remaining common times match the uninterrupted run, but complete
  history equivalence still fails.

`read_parameters` opens those CSVs with `STATUS='unknown'`, without append
positioning. Binary restart skips the initial `output_runout` call, while
version 4 saves masks and schedules but not `x0_runout`, `y0_runout`, initial
runout offsets or previous barycentre history. ETNA's automatic initial origin
therefore remains at its input sentinel after reload; the fissural lava origin
is independently reconstructed from the input, explaining its matching common
values. Previous-barycentre history also participates in steady-stop diagnostics.
Both tested inputs have effective `EPS_STOP=0`; the successful field comparison
must not be extended to arbitrary positive stopping thresholds.

The next correction must preserve the necessary diagnostic history/reference
through restart and preserve or explicitly continue CSV history, with a clear
checkpoint-version compatibility policy. These sources have not been modified
in this package. Passing physical-state continuation alone does not close N9.

## Reproduction and evidence

`TESTS/ACCEPTANCE/n9_long_cases.py` runs the cases outside the repository using
hash-verified executables from an acceptance build manifest. Its completed-run
characterization mode is read-only and retains the failed diagnostic checks;
the command returns nonzero for them. The helper tests verify unique input
overrides, restart record ordering/version/markers and positive/negative CSV
history controls. Commands are documented in the acceptance README.

The committed `results_n9_long_cases_2026_10_11.json` pins source, fixture,
reference and executable checksums, metadata, decoded-field fingerprints,
budget observations and separate physical/diagnostic outcomes. Full runtime
evidence is under `/tmp/imex-n9-long-20261011-v4`; full/split payloads are under
`/tmp/imex-n9-long-20261011-v1`, and preceding lava full runs under
`/tmp/imex-etna-lava-validation-tzq87gsk`. Temporary directories are not a
permanent golden-output archive: reproduce from the pinned inputs/source or
retain those payloads separately before temporary-file cleanup.
