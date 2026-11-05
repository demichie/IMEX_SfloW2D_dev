#!/usr/bin/env python3
"""Enforce the production boundary of the retired state conversion module.

This is a static architectural check, not a numerical equivalence test:
the former state backends have been removed from production.
"""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LEGACY_NAMES = {"state_2d", "state_conversion_2d"}
USE_LEGACY = re.compile(r"^\s*use\s*(?:,\s*(?:non_intrinsic|intrinsic)\s*::\s*|::\s*)?(?:state_2d|state_conversion_2d)\b", re.I)
MAKE_LEGACY = re.compile(r"^\s*spatial/hyperbolic_2d\.\$\(OBJEXT\)\s*:\s*state/state_conversion_2d\.\$\(OBJEXT\)\s*$")
errors = []
for source in sorted(SRC.rglob("*.f90")):
    for lineno, line in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
        if USE_LEGACY.match(line.split("!", 1)[0]):
            errors.append(f"{source.relative_to(ROOT)}:{lineno}: legacy module imported")
for lineno, line in enumerate((SRC / "Makefile.am").read_text(encoding="utf-8").splitlines(), 1):
    if MAKE_LEGACY.match(line):
        errors.append(f"src/Makefile.am:{lineno}: stale legacy prerequisite")
for obsolete in ("state_2d.f90", "state_conversion_2d.f90"):
    if (SRC / "state" / obsolete).exists():
        errors.append(f"src/state/{obsolete}: retired module remains")
if errors:
    print("state_BOUNDARY_FAIL")
    for error in errors:
        print(error)
    sys.exit(1)
print("state_BOUNDARY_PASS: no legacy imports in production; obsolete Make edge absent")
