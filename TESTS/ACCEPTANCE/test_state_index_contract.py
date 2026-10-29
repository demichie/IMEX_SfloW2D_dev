#!/usr/bin/env python3
"""Prevent fixed indices for semantic state components in production Fortran.

Only arrays with a known physical meaning are checked. Compact local vectors,
geometry arrays, and numerical work arrays legitimately use integer positions.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src"
PATTERNS = (
    re.compile(r"\b(?:qp|qpj|qp_source|source_bdry)\s*\(\s*[1-4]\s*[,)]", re.I),
    re.compile(r"\b(?:state%q|simulation%state%q)\s*\(\s*[1-4]\s*,", re.I),
    re.compile(r"\b(?:state%qp|simulation%state%qp|this%qp)\s*\(\s*[1-4]\s*,", re.I),
    re.compile(r"\b(?:eqns_term)\s*\(\s*[1-4]\s*[,)]", re.I),
)
failures: list[str] = []
for source in sorted(SOURCE.rglob("*.f90")):
    if source.name == "kt_ritter_1d.f90":
        continue
    for lineno, original in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
        line = original.split("!", 1)[0]
        for pattern in PATTERNS:
            if pattern.search(line):
                failures.append(f"{source.relative_to(ROOT)}:{lineno}: {original.strip()}")
if failures:
    print("Fixed physical state positions remain:")
    print("\n".join(failures))
    raise SystemExit(1)
print("PASS: semantic primitive and conservative state indices are symbolic")
