# EXAMPLES input compatibility audit — 2026-11-02

The previous D17 static check only checked a small list of deprecated identifiers. The checker now derives the current Fortran NAMELIST vocabulary from `src/io/inpout_2d.f90` and `src/numerics/source_input_2d.f90` and checks assignments in every tracked example template/input against it, including the Etna ensemble input.

**Validation pending:** run `python3 tools/d17/audit_example_inputs.py` on the updated checkout and inspect any reported unknown fields. This static audit does not substitute for generating the inputs and running the solver. No claim is made that all examples have passed on the current revision.

## Confirmed input correction

- `EXAMPLE_DAMBREAK/IMEX_SfloW2D.template`: moved `GAS_LOSS_FLAG=T` out of `PORE_PRESSURE_PARAMETERS` into `MASS_EXCHANGE_PARAMETERS`. The latter is the only current NAMELIST declaring that parameter. Its logical value remains unchanged; the pore-pressure physical coefficients remain unchanged. Generated inputs inherit this fix from the template.

This correction is source-level only until the generated example is parsed and executed with the current solver.
