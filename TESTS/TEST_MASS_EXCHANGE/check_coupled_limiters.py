#!/usr/bin/env python3
"""Independent algebraic regression for the single-carrier exchange budgets.

This tests the limiter mathematics, not the Fortran solver integration.
Neither mixed gas/liquid states nor a separately transported liquid exist.
"""
from itertools import product

PACKING = 0.60
TOL = 1.0e-10


def erosion_factor(depth, alpha, erosion_solid, erosion_carrier, dt):
    deficit = (1.0 - PACKING) * erosion_solid - PACKING * erosion_carrier
    if deficit <= 0.0:
        return 1.0
    return min(1.0, max(0.0, depth * (PACKING - alpha)) / (dt * deficit))


def packing_factor(depth, alpha, erosion_solid, erosion_carrier,
                   deposition_solid, deposition_carrier, carrier_loss, dt):
    margin = PACKING * (depth + dt * (erosion_solid + erosion_carrier))
    margin -= depth * alpha + dt * erosion_solid
    coefficient = PACKING * (deposition_solid + deposition_carrier + carrier_loss)
    coefficient -= deposition_solid
    if coefficient <= 0.0:
        return 1.0
    return min(1.0, max(0.0, margin) / (dt * coefficient))


def check_case(alpha, erosion, deposition, carrier_loss, porosity, carrier):
    depth, dt = 1.0, 1.0
    erosion_carrier = erosion * porosity
    deposition_carrier = 0.10 * deposition

    ef = erosion_factor(depth, alpha, erosion, erosion_carrier, dt)
    erosion *= ef
    erosion_carrier *= ef

    inventory = depth * max(0.0, 1.0 - alpha)
    withdrawal = deposition_carrier + carrier_loss
    factor = min(1.0, inventory / withdrawal) if withdrawal else 1.0
    factor = min(factor, packing_factor(depth, alpha, erosion, erosion_carrier,
                                       deposition, deposition_carrier, carrier_loss, dt))
    deposition *= factor
    deposition_carrier *= factor
    carrier_loss *= factor

    final_solids = depth * alpha + dt * (erosion - deposition)
    final_depth = depth + dt * (erosion + erosion_carrier
                               - deposition - deposition_carrier - carrier_loss)
    assert final_depth > 0.0, carrier
    assert final_solids <= PACKING * final_depth + TOL, carrier
    assert deposition_carrier + carrier_loss <= inventory / dt + TOL, carrier


def main():
    count = 0
    for carrier, alpha, erosion, deposition, loss, porosity in product(
            ("gas", "liquid"), (0.05, 0.25, 0.45, 0.55, 0.60),
            (0.0, 0.005, 0.03, 0.30, 1.0), (0.0, 0.01, 0.10),
            (0.0, 0.01, 0.20), (0.0, 0.10, 1.0)):
        check_case(alpha, erosion, deposition, loss, porosity, carrier)
        count += 1
    print(f"PASS: {count} single-carrier algebraic limiter cases")


if __name__ == "__main__":
    main()
