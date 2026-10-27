#!/usr/bin/env python3
"""D17 cross-example static contract gate. Does not claim solver equivalence."""
import argparse
import ast
import json
from pathlib import Path
import re
import sys

CASES = [
    "EXAMPLE_2D",
    "EXAMPLE_2D_PP",
    "EXAMPLE_BOX",
    "EXAMPLE_BUMP/SUBCRITICAL",
    "EXAMPLE_BUMP/SUPERCRITICAL",
    "EXAMPLE_BUMP/TRANSCRITICAL_WITHOUT_SHOCK",
    "EXAMPLE_BUMP/TRANSCRITICAL_WITH_SHOCK",
    "EXAMPLE_BW",
    "EXAMPLE_BW2D",
    "EXAMPLE_CF",
    "EXAMPLE_DAMBREAK",
    "EXAMPLE_ENGWELL2016",
    "EXAMPLE_ETNA",
    "EXAMPLE_ETNA_LAVA",
    "EXAMPLE_INCLINED_EXCAVATION_2D",
    "EXAMPLE_LATERAL_SOURCE",
    "EXAMPLE_SLOPE",
    "EXAMPLE_TAAL",
    "EXAMPLE_TAAL_BOTTOM"
]
LEGACY_INTERNAL = re.compile(r"\b(?:rheology_model_id|rheology_id_[A-Za-z_0-9]+|source_solver_type|term_ids|SOURCE_TERM_[A-Z_0-9]+|IMEX_VELOCITY_GROUP|IMEX_ENERGY_GROUP|IMEX_PORE_PRESSURE_GROUP|IMEX_GROUP_METHOD|SOURCE_GROUP_ORDER)\b|(?im)^\s*[&$]SOURCE_GROUPS\b", re.IGNORECASE)
NAMELIST = re.compile(r"(?m)^\s*&[A-Za-z_][A-Za-z0-9_]*\s*$")
PLACEHOLDER = re.compile(r"(?<![A-Za-z0-9_])(?:runname|restartfile|nx_cells|ny_cells|x_min|y_min|x_left|dx|n_cells)\b")

def examine(root, name):
    directory = root / "EXAMPLES" / name
    issues = []
    evidence = {"case": name, "checks": {}, "input_type": None, "status": "PASS_STATIC"}
    readme = directory / "README.md"
    if not readme.is_file():
        issues.append("README.md missing")
    else:
        evidence["checks"]["readme"] = "PASS"
    templates = list(directory.glob("IMEX_SfloW2D.template"))
    inputs = list(directory.glob("IMEX_SfloW2D.inp"))
    if not templates and not inputs:
        issues.append("missing namelist template and tracked input")
    for path in templates + inputs:
        text = path.read_text(encoding="utf-8")
        if not NAMELIST.search(text):
            issues.append(f"{path.name}: missing namelist group")
        if LEGACY_INTERNAL.search(text):
            issues.append(f"{path.name}: retired internal identifiers")
        if path.suffix == ".inp" and PLACEHOLDER.search(text):
            issues.append(f"{path.name}: unresolved generator placeholders")
        evidence["checks"][path.name] = "PASS" if not issues else "CHECK"
    evidence["input_type"] = "generated" if templates else "tracked"
    for path in sorted(directory.glob("*.py")):
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (SyntaxError, UnicodeError) as exc:
            issues.append(f"Python syntax: {path.name}: {exc}")
    evidence["checks"]["python_syntax"] = "PASS" if not any("Python syntax:" in x for x in issues) else "FAIL"
    generator = directory / ("create_dambreak.py" if name == "EXAMPLE_DAMBREAK" else "create_example.py")
    if templates and not generator.exists():
        issues.append(f"template present but generator missing: {generator.name}")
    evidence["checks"]["physical_configuration_only"] = "PASS" if not any("retired internal identifiers" in issue for issue in issues) else "FAIL"
    evidence["checks"]["solver_run"] = "NOT_RUN"
    evidence["checks"]["generated_input"] = "NOT_RUN" if templates else "NOT_APPLICABLE"
    if issues:
        evidence["status"] = "FAIL_STATIC"
    evidence["issues"] = issues
    return evidence

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=CASES)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    results = [examine(root, case) for case in (args.case and [args.case] or CASES)]
    result = {"revision": "current checkout", "gate": "static-only", "results": results,
              "summary": {"total": len(results), "pass_static": sum(x["status"] == "PASS_STATIC" for x in results),
                          "fail_static": sum(x["status"] != "PASS_STATIC" for x in results)}}
    payload = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return int(result["summary"]["fail_static"] != 0)

if __name__ == "__main__":
    sys.exit(main())
