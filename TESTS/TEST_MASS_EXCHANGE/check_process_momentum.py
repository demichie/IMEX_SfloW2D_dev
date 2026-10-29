#!/usr/bin/env python3
"""Independent algebraic oracle for the D18.MX per-process momentum contract.

The source of mixture mass includes incoming and outgoing contributions.
Incoming mass carries no initial horizontal momentum; outgoing mass leaves at
the local flow velocity. In particular, zero NET mass does not imply zero
momentum source. This is a reference check, not a Fortran integration test.
"""
from itertools import product
from math import isclose


def source(process_rates, u, v):
    mass = sum(process_rates)
    sink_mass = sum(min(0.0, rate) for rate in process_rates)
    return (mass, u * sink_mass, v * sink_mass)


def expected_from_independent_budgets(incoming, outgoing, u, v):
    """Outgoing is a positive, physically removed inventory rate."""
    return (incoming - outgoing, -outgoing * u, -outgoing * v)


count = 0
for incoming, outgoing, u, v in product(
    (0.0, 0.5, 2.0, 10.0),
    (0.0, 0.5, 2.0, 10.0),
    (-3.0, 0.0, 1.5),
    (-2.0, 0.0, 4.0),
):
    # Separate representative incoming erosion/entrainment and outgoing
    # deposition/loss; partitioning does not alter the final momentum.
    rates = (0.25 * incoming, 0.75 * incoming,
             -0.4 * outgoing, -0.6 * outgoing)
    actual = source(rates, u, v)
    expected = expected_from_independent_budgets(incoming, outgoing, u, v)
    for observed, reference in zip(actual, expected):
        assert isclose(observed, reference, rel_tol=1e-13, abs_tol=1e-13), (
            rates, u, v, actual, expected
        )
    if incoming == outgoing and outgoing > 0.0 and (u != 0.0 or v != 0.0):
        assert actual[0] == 0.0
        assert actual[1] != 0.0 or actual[2] != 0.0
    count += 1

print(f"PASS: D18.MX independent per-process momentum oracle ({count} cases)")
