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
    for failure in failures:
        print("FAIL:", failure, file=sys.stderr)
    print(f"D17 physical treatment audit: {len(catalog)} configured sources, "
          f"{len(catalog) * len(METHODS)} concrete interfaces; "
          f"{len(failures)} failures")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
