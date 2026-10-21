"""Check the N8 fixture's production-update geometry diagnostic records."""

import math
import hashlib
import json
from pathlib import Path
import sys


def check(path, profile=None, requested=None):
    actual_threads = None
    evidence = []
    flat_cases = []
    for line in path.read_text().splitlines():
        if "Actual OpenMP team:" in line:
            actual_threads = int(line.rsplit(":", 1)[1])
        elif "N8_FLAT_EQUIVALENCE " in line:
            flat_cases.append(line.split()[-1])
        elif line.startswith("N8_EXPECTED_DIAGNOSTICS "):
            fields = line.split()
            if len(fields) != 7:
                raise AssertionError("incorrect diagnostic field count")
            values = [float(value) for value in fields[2:]]
            if not all(math.isfinite(value) for value in values):
                raise AssertionError("nonfinite production diagnostic")
            cell_volume, geometric_volume, mismatch, mismatch_l1, mismatch_linf = values
            scale = max(1.0, abs(cell_volume), abs(geometric_volume), abs(mismatch))
            if (abs(geometric_volume-cell_volume-mismatch) > 2e-12*scale
                    or abs(mismatch) > mismatch_l1 + 2e-12*scale
                    or mismatch_l1 < 0.0 or mismatch_linf < 0.0):
                raise AssertionError("inconsistent production geometry diagnostics")
            evidence.append({"case": fields[1], "diagnostics": values})
    if not evidence:
        raise AssertionError("no production diagnostic records")
    print(f"PASS: {len(evidence)} production bed-volume/mismatch diagnostic records")
    if profile is not None:
        required = ["deposition", "erosion", "combined", "masked", "masked_fissural",
                    "cutoff", "loss_only", "loss_reserve", "thermal", "restart_first",
                    "restart_continuous", "restart_resumed", "gas_packing", "gas_reserve",
                    "gas_exhausted", "gas_inflow", "gas_combined", "gas_masked"]
        if ([item["case"] for item in evidence] != required or flat_cases != required
                or actual_threads != requested):
            raise AssertionError("incomplete N8 case set or incorrect actual OpenMP team")
        hashes = {}
        for filename in ("snapshot.bin", "checkpoint.bin", "gas_snapshot.bin", "flat_evaluations.bin"):
            hashes[filename] = hashlib.sha256((path.parent / filename).read_bytes()).hexdigest()
        print("N8_PRODUCTION_EVIDENCE " + json.dumps({
            "profile": profile, "actual_threads": actual_threads,
            "sha256": hashes, "diagnostics": evidence, "flat_cases": flat_cases,
        }, sort_keys=True))


if __name__ == "__main__":
    check(Path(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else None,
          int(sys.argv[3]) if len(sys.argv) > 3 else None)
