"""Check production bed-rate diagnostics, including area/volume units."""

import math
import hashlib
import json
from pathlib import Path
import sys


def check(path, profile=None, requested=None):
    volume = mismatch = None
    checks = 0
    actual_threads = None
    evidence = []
    for line in path.read_text().splitlines():
        if "Actual OpenMP team:" in line:
            actual_threads = int(line.rsplit(":", 1)[1])
        elif "evolving-topography volume-rate cell/geometric:" in line:
            volume = [float(v) for v in line.split(":", 1)[1].split()]
        elif "evolving-topography mismatch integral/L1/Linf:" in line:
            mismatch = [float(v) for v in line.split(":", 1)[1].split()]
        elif line.startswith("N8_EXPECTED_DIAGNOSTICS "):
            expected = [float(v) for v in line.split()[2:]]
            if volume is None or mismatch is None:
                raise AssertionError("missing production diagnostics")
            actual = volume + mismatch
            if len(actual) != 5 or len(expected) != 5:
                raise AssertionError("incorrect diagnostic field count")
            for value, reference in zip(actual, expected):
                if not math.isfinite(value) or abs(value-reference) > 2e-12*max(1.0, abs(reference)):
                    raise AssertionError(f"diagnostic {value} != expected {reference}")
            checks += 1
            evidence.append({"case": line.split()[1], "production": actual, "reference": expected})
            volume = mismatch = None
    if not checks:
        raise AssertionError("no diagnostic comparisons executed")
    print(f"PASS: {checks} production bed-volume/mismatch diagnostic comparisons")
    if profile is not None:
        if checks != 12 or actual_threads != requested:
            raise AssertionError("incomplete N8 case set or incorrect actual OpenMP team")
        hashes = {}
        for filename in ("snapshot.bin", "checkpoint.bin"):
            hashes[filename] = hashlib.sha256((path.parent / filename).read_bytes()).hexdigest()
        print("N8_PRODUCTION_EVIDENCE " + json.dumps({
            "profile": profile, "actual_threads": actual_threads,
            "sha256": hashes, "diagnostics": evidence,
        }, sort_keys=True))


if __name__ == "__main__":
    check(Path(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else None,
          int(sys.argv[3]) if len(sys.argv) > 3 else None)
