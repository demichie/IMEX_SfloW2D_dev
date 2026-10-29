#!/usr/bin/env python3
"""Check canonical state-index migration invariants without rewriting code."""
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parents[2]
required = {
    "TESTS/TEST_STATE_CONVERSION/test_state_conversion.f90": "CALL configure_canonical_state_indices()",
    "TESTS/TEST_RHEOLOGY_FROUDE/test_rheology_froude.f90": "CALL configure_canonical_state_indices()",
    "TESTS/TEST_SPATIAL_OPERATOR/test_spatial_operator.f90": "CALL configure_canonical_state_indices()",
    "TESTS/TEST_IMEX_STAGES/test_imex_stages.f90": "CALL configure_canonical_state_indices()",
    "TESTS/TEST_MASS_EXCHANGE/test_mass_exchange.f90": "CALL configure_canonical_state_indices()",
    "TESTS/TEST_MASS_EXCHANGE/test_projection_stress.f90": "CALL configure_canonical_state_indices()",
    "TESTS/ACCEPTANCE/test_d17_state_numeric.f90": "CALL configure_canonical_state_indices()",
}
issues = []
for name, marker in required.items():
    content = (root / name).read_text()
    if marker not in content:
        issues.append(f"{name}: missing canonical fixture configuration")

state = (root / "src/state")
for name in ("state_conversion_2d.f90", "state_conversion_phys_var.inc",
             "velocity_from_conservative.inc"):
    path = state / name
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if line.lstrip().startswith("!"):
            continue
        if re.search(r"\b(?:qc|qj)\s*\(\s*[1-4]\s*\)", line, re.I):
            issues.append(f"{path.relative_to(root)}:{number}: hard-coded conservative index")
if issues:
    print("\n".join(issues), file=sys.stderr)
    sys.exit(1)
print("PASS: canonical state index migration checks")
