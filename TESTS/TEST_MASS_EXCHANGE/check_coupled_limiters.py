#!/usr/bin/env python3
"""Independent algebraic matrix for mass-exchange carrier and packing limiters.

Run with python3 TESTS/TEST_MASS_EXCHANGE/check_coupled_limiters.py.
This checks finite-volume balances, not the Fortran solver integration.
"""
from itertools import product

PACKING = 0.60
TOL = 1.0e-10


def packing_factor(depth, solid_fraction, erosion_solid, erosion_carrier,
                   deposition_solid, deposition_carrier, shared_loss,
                   other_loss, dt):
    """Independent solution of final-solid <= packing * final-depth."""
    margin = PACKING * (depth + dt * (erosion_solid + erosion_carrier - other_loss))
    margin -= depth * solid_fraction + dt * erosion_solid
    coefficient = PACKING * (deposition_solid + deposition_carrier + shared_loss)
    coefficient -= deposition_solid
    if coefficient <= 0:
        return 1.0
    return min(1.0, max(0.0, margin) / (dt * coefficient))


def check_case(gas, liquid, alpha, alpha_liquid, erosion, deposition,
               gas_loss, liquid_loss):
    """Check each admissible carrier configuration with simultaneous processes."""
    depth = 1.0
    dt = 1.0
    erosion_carrier = 0.10 * erosion
    deposition_carrier = 0.10 * deposition

    if gas and liquid:
        liquid_scale = min(1.0, alpha_liquid / liquid_loss) if liquid_loss else 1.0
        liquid_scale = min(liquid_scale, packing_factor(
            depth, alpha, erosion, erosion_carrier, 0, 0, liquid_loss, 0, dt))
        liquid_loss *= liquid_scale

    if gas:
        gas_inventory = max(0.0, depth * (1.0 - alpha - alpha_liquid))
        total_gas = deposition_carrier + gas_loss
        gas_scale = min(1.0, gas_inventory / total_gas) if total_gas else 1.0
        gas_scale = min(gas_scale, packing_factor(
            depth, alpha, erosion, erosion_carrier, deposition,
            deposition_carrier, gas_loss, liquid_loss, dt))
        deposition *= gas_scale
        deposition_carrier *= gas_scale
        gas_loss *= gas_scale
    else:
        total_liquid = deposition_carrier + liquid_loss
        liquid_scale = min(1.0, (1.0 - alpha) / total_liquid) if total_liquid else 1.0
        liquid_scale = min(liquid_scale, packing_factor(
            depth, alpha, erosion, erosion_carrier, deposition,
            deposition_carrier, liquid_loss, 0, dt))
        deposition *= liquid_scale
        deposition_carrier *= liquid_scale
        liquid_loss *= liquid_scale

    final_solids = alpha + dt * (erosion - deposition)
    final_depth = depth + dt * (
        erosion + erosion_carrier - deposition - deposition_carrier
        - gas_loss - liquid_loss)
    assert final_depth > 0
    assert final_solids <= PACKING * final_depth + TOL
    if gas:
        assert gas_loss + deposition_carrier <= (1.0 - alpha - alpha_liquid) + TOL
    if gas and liquid:
        assert liquid_loss <= alpha_liquid + TOL
    if not gas:
        assert deposition_carrier + liquid_loss <= (1.0 - alpha) + TOL


def main():
    count = 0
    for gas, liquid in [(True, False), (True, True), (False, True)]:
        for alpha, erosion, deposition, gas_loss, liquid_loss in product(
                [0.05, 0.25, 0.45, 0.55], [0.0, 0.005, 0.03],
                [0.0, 0.01, 0.1], [0.0, 0.01, 0.2], [0.0, 0.01, 0.2]):
            if not gas:
                gas_loss = 0.0
            if not liquid:
                liquid_loss = 0.0
            alpha_liquid = 0.2 if gas and liquid else 0.0
            check_case(gas, liquid, alpha, alpha_liquid, erosion,
                       deposition, gas_loss, liquid_loss)
            count += 1
    print(f"PASS: {count} independent algebraic mass-exchange limiter cases")


if __name__ == "__main__":
    main()
