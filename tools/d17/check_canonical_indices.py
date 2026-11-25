#!/usr/bin/env python3
"""Check symbolic state indices in the current canonical state modules."""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
REQUIRED_FIXTURES = (
    "TESTS/TEST_STATE/test_state.f90",
    "TESTS/TEST_STATE_CONVERSION/test_state_conversion.f90",
    "TESTS/TEST_RHEOLOGY_FROUDE/test_rheology_froude.f90",
    "TESTS/TEST_SPATIAL_OPERATOR/test_spatial_operator.f90",
    "TESTS/TEST_MASS_EXCHANGE/test_projection_stress.f90",
)
FIXTURE_CALL = "CALL configure_canonical_state_indices()"
FIXED_INDEX = re.compile(r"\b(?:qc|qj)\s*\(\s*[1-4]\s*\)", re.I)


def main() -> int:
    issues = []
    for name in REQUIRED_FIXTURES:
        path = ROOT / name
        if not path.is_file():
            issues.append(f"{name}: test fixture missing")
        elif FIXTURE_CALL not in path.read_text(encoding="utf-8"):
            issues.append(f"{name}: missing canonical fixture configuration")

    state_dir = ROOT / "src/state"
    state_sources = sorted(state_dir.glob("*.f90"))
    state_includes = sorted((state_dir / "inc").glob("*.inc"))
    if not state_sources or not state_includes:
        issues.append("canonical state modules or shared kernels missing")
    for path in state_sources + state_includes:
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if line.lstrip().startswith("!"):
                continue
            if FIXED_INDEX.search(line.split("!", 1)[0]):
                issues.append(f"{path.relative_to(ROOT)}:{number}: hard-coded conservative index")
    if issues:
        print("\n".join(issues), file=sys.stderr)
        return 1
    print(f"PASS: canonical indices in {len(state_sources)} state modules, "
          f"{len(state_includes)} shared kernels and {len(REQUIRED_FIXTURES)} fixtures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
