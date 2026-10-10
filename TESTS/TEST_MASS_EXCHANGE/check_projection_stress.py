"""Frozen N8-B case inventory, refinement limits and separately scoped area algebra."""
import hashlib
import json
import math
from pathlib import Path
import sys

SIZES = ((1, 16, 1), (1, 32, 1), (1, 64, 1),
         (2, 16, 12), (2, 32, 24), (2, 64, 48))
STEPS = {"smooth": 1, "wet_dry": 1, "rough": 32, "repeated": 64}
ROUND_OFF = 2048 * sys.float_info.epsilon
MAX_REFINEMENT_RATIO = 0.8


def unequal_area_algebra():
    """Check AC/4 mass-lumping algebra, not nonuniform production support."""
    areas = ((1.0, 2.0, 3.0), (0.75, 4.0, 1.5))
    fields = {"constant": ((2.0, 2.0, 2.0), (2.0, 2.0, 2.0)),
              "signed": ((1.0, -2.0, 3.0), (-4.0, 5.0, -6.0)),
              "isolated": ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))}
    records = []
    wrong_weights_rejected = False
    for name, field in fields.items():
        weights = [[0.0]*4 for _ in range(3)]
        loads = [[0.0]*4 for _ in range(3)]
        for k in range(2):
            for j in range(3):
                for b in (k, k+1):
                    for a in (j, j+1):
                        weights[b][a] += areas[k][j]/4
                        loads[b][a] += areas[k][j]*field[k][j]/4
        nodes = [[loads[k][j]/weights[k][j] for j in range(4)] for k in range(3)]
        cell_integral = sum(areas[k][j]*field[k][j] for k in range(2) for j in range(3))
        node_integral = sum(weights[k][j]*nodes[k][j] for k in range(3) for j in range(4))
        # Independently integrate the Q1 nodal field on each unequal-area cell.
        q1_integral = sum(areas[k][j]*sum(nodes[b][a] for b in (k, k+1)
                           for a in (j, j+1))/4 for k in range(2) for j in range(3))
        scale = max(1.0, sum(abs(areas[k][j]*field[k][j]) for k in range(2) for j in range(3)))
        if max(abs(node_integral-cell_integral), abs(q1_integral-cell_integral)) > ROUND_OFF*scale:
            raise AssertionError("unequal-area mass-lumped algebra")
        arithmetic_nodes = [[sum(field[b][a] for b in range(max(0, k-1), min(2, k+1))
                             for a in range(max(0, j-1), min(3, j+1))) /
                            ((min(2, k+1)-max(0, k-1))*(min(3, j+1)-max(0, j-1)))
                            for j in range(4)] for k in range(3)]
        wrong_integral = sum(weights[k][j]*arithmetic_nodes[k][j] for k in range(3) for j in range(4))
        wrong_weights_rejected |= abs(wrong_integral-cell_integral) > ROUND_OFF*scale
        records.append({"case": name, "cell_integral": cell_integral,
                        "node_integral": node_integral, "q1_integral": q1_integral,
                        "unweighted_integral": wrong_integral})
    if not wrong_weights_rejected:
        raise AssertionError("unequal-area negative control not exercised")
    return {"scope": "synthetic algebra only; production remains uniform Cartesian",
            "areas": areas, "records": records, "wrong_weights_rejected": True}


def check(path, profile, requested):
    records = []
    actual = None
    for line in path.read_text().splitlines():
        if "Actual OpenMP team:" in line:
            actual = int(line.rsplit(":", 1)[1])
        elif line.startswith("N8B_RECORD "):
            fields = line.split()
            if len(fields) != 12:
                raise AssertionError("incorrect stress record size")
            dimension, nx, ny, iteration = map(int, fields[2:6])
            metrics = list(map(float, fields[6:]))
            if not all(math.isfinite(v) for v in metrics):
                raise AssertionError("nonfinite stress diagnostic")
            volume_cell, volume_geom, signed, l1, linf, outside = metrics
            if max(abs(volume_cell-volume_geom), abs(signed), abs(outside)) > ROUND_OFF:
                raise AssertionError("volume identity or front localization")
            if min(l1, linf) <= ROUND_OFF:
                raise AssertionError("stress not exercised")
            records.append({"case": fields[1], "dimension": dimension,
                            "nx": nx, "ny": ny, "iteration": iteration, "metrics": metrics})
    expected = [(name, d, nx, ny, i) for d, nx, ny in SIZES
                for name, steps in STEPS.items() for i in range(1, steps+1)]
    observed = [(r["case"], r["dimension"], r["nx"], r["ny"], r["iteration"]) for r in records]
    if observed != expected or actual != requested:
        raise AssertionError("incomplete stress case set or incorrect actual team")
    refinements = []
    for dimension in (1, 2):
        for name in ("smooth", "wet_dry"):
            series = [r for r in records if r["dimension"] == dimension and r["case"] == name]
            ratios = [b["metrics"][3]/a["metrics"][3] for a, b in zip(series, series[1:])]
            if not all(0 < ratio < MAX_REFINEMENT_RATIO for ratio in ratios):
                raise AssertionError(f"{name} dimension {dimension}: L1 does not decrease")
            refinements.append({"dimension": dimension, "case": name, "L1_ratios": ratios,
                                "L1": [r["metrics"][3] for r in series],
                                "Linf": [r["metrics"][4] for r in series]})
    evidence = {"profile": profile, "actual_threads": actual, "records": records,
                "refinements": refinements, "unequal_area_algebra": unequal_area_algebra(),
                "snapshot_sha256": hashlib.sha256((path.parent/"projection_snapshot.bin").read_bytes()).hexdigest(),
                "log_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    print("N8B_PRODUCTION_EVIDENCE " + json.dumps(evidence, sort_keys=True))
    return evidence


if __name__ == "__main__":
    check(Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]))
