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


def structural_warnings(path):
    """Conservative nesting diagnostics for free-form IF/DO/SELECT/WHERE blocks.

    Only inspect full (non-continued) control statements. Leave generated
    directives and multiline conditions to the compiler and human review.
    """
    import re
    stack = []
    warnings = []
    previous_continued = False
    for line_no, line in enumerate(path.read_text().splitlines(), 1):
        text = line.strip()
        if not text or text.startswith(("!", "#")):
            continue
        if previous_continued or text.startswith("&"):
            previous_continued = text.endswith("&")
            continue
        previous_continued = text.endswith("&")
        if previous_continued:
            continue
        prefix = len(line) - len(line.lstrip(" "))
        upper = text.upper()
        opener = None
        if re.match(r"^IF\s*\(.*\)\s*THEN\s*$", upper):
            opener = "IF"
        elif re.match(r"^DO(?:\s|$)", upper):
            opener = "DO"
        elif re.match(r"^SELECT\s+(?:CASE|TYPE|RANK)\s*\(", upper):
            opener = "SELECT"
        elif re.match(r"^WHERE\s*\(.*\)\s*$", upper):
            opener = "WHERE"
        closer = re.match(r"^END\s*(IF|DO|SELECT|WHERE)\b", upper)
        if closer:
            kind = closer.group(1)
            if stack and stack[-1][0] == kind:
                start_kind, start_indent, start_line = stack.pop()
                if prefix != start_indent:
                    warnings.append(
                        f"{path}:{line_no}: END {kind} indent {prefix}, "
                        f"opener line {start_line} indent {start_indent}")
        elif re.match(r"^(ELSE(?:\s*IF)?\b|CASE\b)", upper):
            kind = "IF" if upper.startswith("ELSE") else "SELECT"
            if stack and stack[-1][0] == kind and prefix != stack[-1][1]:
                warnings.append(
                    f"{path}:{line_no}: branch indent {prefix}, "
                    f"opening {kind} indent {stack[-1][1]}")
        if opener:
            if stack and prefix <= stack[-1][1]:
                warnings.append(
                    f"{path}:{line_no}: nested {opener} should be indented "
                    f"more than enclosing {stack[-1][0]} at line {stack[-1][2]}")
            stack.append((opener, prefix, line_no))
    return warnings

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="Fortran files/directories")
    parser.add_argument("--structural-report", action="store_true",
                        help="print conservative block-nesting warnings (non-failing)")
    args = parser.parse_args()
    errors = []
    seen = set()
    nesting_warnings = []
    for path in files(args.paths):
        if path in seen:
            continue
        seen.add(path)
        if args.structural_report:
            nesting_warnings.extend(structural_warnings(path))
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
    if args.structural_report:
        print(*nesting_warnings, sep=chr(10)) if nesting_warnings else print('No conservative nesting warnings')
        print(f'Structural warnings: {len(nesting_warnings)}')
    if errors:
        print("\n".join(errors))
        return 1
    print(f"PASS: Fortran whitespace and four-space indentation in {len(seen)} files")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
