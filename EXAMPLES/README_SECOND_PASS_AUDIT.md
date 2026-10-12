# EXAMPLES README migration — second-pass assessment (2026-10-12)

Source repository: `demichie/IMEX_SfloW2D_dev` at `2bb0efff87ee8f8fb93815fc7e050011c3949acf`.

This archive is extracted **inside** `EXAMPLES/` and contains 20 standardized `README.md` files. It deliberately does **not** delete or replace any of the original `README.txt` files, scripts, input namelists or DEM data. The legacy text retains the archival wording and any fine detail missing from a standardized summary; relative links to each legacy README were added explicitly.

## Validation carried out

- Reviewed the 20 original README texts through the GitHub repository interface, including nested BUMP variants.
- Compared the generated documentation with the legacy descriptions; inspected the generator arguments for the relevant complex cases and the committed runtime inputs for ETNA_LAVA, TAAL, TAAL_BOTTOM and CF.
- Checked the YAML metadata, required headings and balanced Markdown code fences for all 20 standardized READMEs.
- Preserved bibliographic references in the first-pass README sections and retained access to all original documentation.

## Corrections and caveats

- `EXAMPLE_BW2D`: eight script arguments; conflicting pulse narratives not silently erased.
- `EXAMPLE_LATERAL_SOURCE`: source mass-flow rate, not initial thickness; units remain to be confirmed from the entire generator.
- `EXAMPLE_ENGWELL2016`: script uses `xs_init`; inspect transformation before certifying physical interpretation.
- `EXAMPLE_DAMBREAK`: current generator `create_dambreak.py`, with older filenames identified as historical; authoritative descriptions of pore pressure and inhibited gas escape retained.
- `EXAMPLE_ETNA`: `cleanFolder.sh` capitalization and ellipsoid-preprocessing explanation.
- `EXAMPLE_2D`: old `create_example4.py` is not the committed generator filename.
- `EXAMPLE_ENSEMBLE_ETNA`: placeholder `ENSEMBLE_`, link/template execution stages retained.
- `EXAMPLE_CF`, `EXAMPLE_TAAL` and `EXAMPLE_TAAL_BOTTOM`: their committed namelists have ESRI output disabled, contrary to the generic legacy output descriptions.
- `EXAMPLE_ETNA_LAVA`: extra DEM and acceptance context restored; validation relates to a pinned source and input configuration.

## Outstanding work before deleting README.txt

This pass **does not certify lossless transfer of every original sentence**. The archived `README.txt` files must remain until a byte/section-level completeness comparison of all paragraphs, parameter tables, references and usage commands is completed in a local working checkout. Furthermore, all 20 runs have **not** been executed from freshly generated fixtures. Some generated-input parameters, postprocessing dependencies and source mass/volume interpretation are not yet fully verified. Do not claim execution validation based only on README inspection.

Recommended final step: run a strict local comparison script against the original `README.txt`, manually sign off every change, and only then optionally remove the historical files. The MkDocs generator should render the entire Markdown body rather than extracting only required sections.
