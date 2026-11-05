#!/usr/bin/env python3
"""Enforce the production boundary of the retired state conversion module.

This is a static architectural check, not a numerical equivalence test:
the legacy module remains built for the old/new characterization tests.
"""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LEGACY = SRC / "state" / "state_conversion_2d.f90"
USE_LEGACY = re.compile(r"^\s*use\s*(?:,\s*(?:non_intrinsic|intrinsic)\s*::\s*|::\s*)?state_conversion_2d\b", re.I)
MAKE_LEGACY = re.compile(r"^\s*spatial/hyperbolic_2d\.\$\(OBJEXT\)\s*:\s*state/state_conversion_2d\.\$\(OBJEXT\)\s*$")
errors = []
for source in sorted(SRC.rglob("*.f90")):
    if source == LEGACY:
        continue
    for lineno, line in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
        if USE_LEGACY.match(line.split("!", 1)[0]):
            errors.append(f"{source.relative_to(ROOT)}:{lineno}: legacy module imported")
for lineno, line in enumerate((SRC / "Makefile.am").read_text(encoding="utf-8").splitlines(), 1):
    if MAKE_LEGACY.match(line):
        errors.append(f"src/Makefile.am:{lineno}: stale legacy prerequisite")
if errors:
    print("STATE_NEW_BOUNDARY_FAIL")
    for error in errors:
        print(error)
    sys.exit(1)
print("STATE_NEW_BOUNDARY_PASS: no legacy imports in production; obsolete Make edge absent")
