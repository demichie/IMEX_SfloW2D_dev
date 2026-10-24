#!/usr/bin/env python3
"""Validate concrete four-treatment entry points for every D17 catalog source.

A source treatment must have a real Fortran procedure with a consistent name
in src/physics/. A public generic alias is not sufficient.
"""
from __future__ import annotations

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "src/numerics/source_catalog_2d.f90"
PHYSICS = ROOT / "src/physics"
METHODS = ("explicit", "semi_implicit", "implicit", "analytic")


def procedure_body(code: str, name: str) -> str | None:
    escaped = re.escape(name)
    header = re.search(
        rf"(?im)^\s*(?:pure\s+)?(?:real\(wp\)\s+|complex\(wp\)\s+)?"
        rf"(?:subroutine|function)\s+{escaped}\s*\(",
        code,
    )
    if not header:
        return None
    end = re.search(
        rf"(?im)^\s*end\s+(?:subroutine|function)\s+{escaped}\b",
        code[header.end():],
    )
    if not end:
        return None
    return code[header.end():header.end() + end.start()]


def main() -> int:
    text = CATALOG.read_text(encoding="utf-8")
    catalog = re.findall(
        r"source_term_type\('([^']+)'\s*,\s*'([^']+)'\s*,"
        r"\s*\.(true|false)\.\s*,\s*\.(true|false)\.\s*,"
        r"\s*\.(true|false)\.\s*\)",
        text,
        flags=re.IGNORECASE,
    )
    if not catalog:
        print("ERROR: no D17 catalog entries found", file=sys.stderr)
        return 1

    code_by_file = {
        f: f.read_text(encoding="utf-8")
        for f in PHYSICS.rglob("*.f90")
    }
    failures = []
    for source, target, explicit, newton, analytic in catalog:
        supported = {
            "explicit": explicit.lower() == "true",
            "semi_implicit": newton.lower() == "true",
            "implicit": newton.lower() == "true",
            "analytic": analytic.lower() == "true",
        }
        for method in METHODS:
            proc = f"{'solve' if method == 'analytic' else 'eval'}_{source}_{method}"
            owners = [
                (file, contents, procedure_body(contents, proc))
                for file, contents in code_by_file.items()
            ]
            owners = [(file, contents, body) for file, contents, body in owners if body is not None]
            if len(owners) != 1:
                failures.append(f"{source}: {proc} has {len(owners)} physical implementations (expected 1)")
                continue
            file, contents, body = owners[0]
            public = re.search(rf"(?im)^\s*public\s*::[^\n]*\b{re.escape(proc)}\b", contents)
            if not public:
                failures.append(f"{source}: {proc} is not PUBLIC in {file.relative_to(ROOT)}")
            rejects_method = bool(re.search(r"(?i)\bcall\s+fatal_error\s*\(\s*['\"][^'\"]*unsupported", body))
            if supported[method] and rejects_method:
                failures.append(f"{source}: supported {method} uses fatal_error in {file.relative_to(ROOT)}")
            if not supported[method] and not rejects_method:
                failures.append(f"{source}: unsupported {method} must reject accidental invocation")
    # Each physical source owns a single, ordered four-entry public contract.
    # Keeping all four procedures together also prevents accidentally routing
    # different treatments through unrelated physical modules.
    for source, _target, _explicit, _newton, _analytic in catalog:
        names = [
            f"{'solve' if method == 'analytic' else 'eval'}_{source}_{method}"
            for method in METHODS
        ]
        owners = []
        for name in names:
            matches = [
                (file, re.search(
                    rf"(?im)^\s*(?:pure\s+)?(?:real\(wp\)\s+|complex\(wp\)\s+)?"
                    rf"(?:subroutine|function)\s+{re.escape(name)}\s*\(",
                    contents,
                ))
                for file, contents in code_by_file.items()
            ]
            owners.append([(file, m.start()) for file, m in matches if m])
        if not all(len(matches) == 1 for matches in owners):
            continue  # Existing checks already report missing/duplicated procedures.
        if len({matches[0][0] for matches in owners}) != 1:
            failures.append(f"{source}: treatment procedures are split across physical modules")
            continue
        positions = [matches[0][1] for matches in owners]
        if positions != sorted(positions):
            failures.append(f"{source}: procedure definitions must be ordered "
                            "explicit, semi_implicit, implicit, analytic")
        file = owners[0][0][0]
        contents = code_by_file[file]
        declarations = []
        for name in names:
            matches = list(re.finditer(
                rf"(?im)^\s*public\s*::[^\n]*\b{re.escape(name)}\b",
                contents,
            ))
            declarations.append(matches)
        if all(len(matches) == 1 for matches in declarations):
            positions = [matches[0].start() for matches in declarations]
            if positions != sorted(positions):
                failures.append(f"{source}: PUBLIC declarations must be ordered "
                                "explicit, semi_implicit, implicit, analytic")
    # D17.14/15: STEP scheduling must reach physical rates exclusively through
    # the common source_terms bridge; inventories remain in the STEP executor.
    executor = (ROOT / "src/numerics/step_source_executor_2d.f90").read_text(
        encoding="utf-8"
    )
    bridge = (ROOT / "src/physics/source_terms_2d.f90").read_text(
        encoding="utf-8"
    )
    equation_bridge = (ROOT / "src/physics/equation_terms_2d.f90").read_text(
        encoding="utf-8"
    )
    if not re.search(
        r"(?im)^\s*use\s+source_terms_2d\s*,\s*only\s*:\s*eval_step_source_explicit\b",
        executor,
    ):
        failures.append("STEP executor must import the common explicit source bridge")
    if len(re.findall(r"(?i)\bcall\s+eval_step_source_explicit\s*\(", executor)) != 1:
        failures.append("STEP executor must evaluate the explicit source bridge once per cell")
    if re.search(r"(?im)^\s*use\s+(?:equation_terms_2d|mass_exchange_assembly_2d)\b",
                 executor):
        failures.append("STEP executor must not import the old equation/mass-rate bridge")
    if re.search(r"(?i)\bcall\s+(?:eval_expl_terms|eval_mass_exchange_terms)\s*\(",
                 executor):
        failures.append("STEP executor must not call the old physical rate entry point")
    if not re.search(r"(?i)\bcall\s+eval_mass_exchange_terms\s*\(", bridge):
        failures.append("common STEP bridge must delegate to the atomic mass-exchange assembly")
    if re.search(r"(?i)\bcall\s+eval_mass_exchange_terms\s*\(",
                 equation_bridge):
        failures.append("equation_terms must not call the mass-exchange kernel directly")
    # D17 cleanup gate: callers must use the canonical source assemblers.
    # Check real USE statements (rather than comments) across the Fortran tree.
    removed = {
        "eval_expl_terms", "eval_stage_expl_terms", "eval_step_expl_terms",
        "eval_semi_implicit_terms", "eval_implicit_terms", "eval_analytic_terms",
    }
    for filename in (ROOT / "src").rglob("*.f90"):
        code = filename.read_text(encoding="utf-8")
        no_comments = "\n".join(line.split("!")[0] for line in code.splitlines())
        for old_name in removed:
            if re.search(rf"(?i)\\b(?:call|module\\s+procedure)\\s+{old_name}\\b", no_comments):
                failures.append(f"legacy source assembler {old_name} in {filename.relative_to(ROOT)}")
        if filename != ROOT / "src/physics/equation_terms_2d.f90":
            for match in re.finditer(
                r"(?im)^\\s*use\\s+equation_terms_2d\\s*,\\s*only\\s*:\\s*"
                r"((?:[^\\n]|\\n\\s*&)*?)(?=\\n\\s*(?:use|implicit|integer|real|type|logical|character|contains|end|public|private|!|$))",
                no_comments,
            ):
                exposed = match.group(1).lower()
                if re.search(r"eval_(?:boundary|bottom|radial|geometry)_(?:injection|curvature)", exposed):
                    failures.append(f"legacy source import through equation_terms in {filename.relative_to(ROOT)}")
    equation_file = ROOT / "src/physics/equation_terms_2d.f90"
    equation_code = equation_file.read_text(encoding="utf-8")
    if re.search(r"(?im)^\s*(?:public\s*::|use\s+rheology_2d\s*,\s*only\s*:)\s*curvature_acceleration\b", equation_code):
        failures.append("equation_terms must not re-export curvature_acceleration")

    # STEP rate/update paths must reuse the same physical dry-state initializer
    # as the RK halo/stage path; do not duplicate ambient-temperature defaults.
    step_code = (ROOT / "src/numerics/step_source_executor_2d.f90").read_text(encoding="utf-8")
    if step_code.count("CALL initialize_dry_physical_cell(") != 2:
        failures.append("STEP executor must reuse dry-cell initializer in both branches")
    if re.search(r"(?im)^\s*qp(?:_cell|\([^\n]*\))?\s*\(4[^\n]*=\s*T_ambient", step_code):
        failures.append("STEP executor must not duplicate ambient dry initialization")
    for failure in failures:
        print("FAIL:", failure, file=sys.stderr)
    print(f"D17 physical treatment audit: {len(catalog)} configured sources, "
          f"{len(catalog) * len(METHODS)} concrete interfaces; "
          f"{len(failures)} failures")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
