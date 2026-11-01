#!/usr/bin/env python3
"""Verify mandatory Fortran whitespace rules; pass paths or directories."""
import argparse
from pathlib import Path

def files(targets):
    for target in targets:
        p = Path(target)
        if p.is_dir():
            yield from sorted(p.rglob("*.f90"))
        elif p.suffix.lower() == ".f90":
            yield p

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="Fortran files/directories")
    args = parser.parse_args()
    errors = []
    seen = set()
    for path in files(args.paths):
        if path in seen:
            continue
        seen.add(path)
        for line_no, line in enumerate(path.read_text().splitlines(), 1):
            if "\t" in line:
                errors.append(f"{path}:{line_no}: TAB forbidden")
            if line.rstrip(" \t") != line:
                errors.append(f"{path}:{line_no}: trailing whitespace")
            if not line.strip():
                continue
            prefix = len(line) - len(line.lstrip(" "))
            if not line.lstrip().startswith("#") and prefix % 4:
                errors.append(f"{path}:{line_no}: indentation is not a multiple of 4")
    if errors:
        print("\n".join(errors))
        return 1
    print(f"PASS: Fortran whitespace and four-space indentation in {len(seen)} files")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
