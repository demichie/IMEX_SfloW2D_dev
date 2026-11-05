#!/usr/bin/env python3
"""Static Fortran module, source-API, and Automake dependency audit.

No compiler needed. Paths are relative to repository root. The audit checks
direct USE edges, not merely an achievable transitive build order.
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src"
MAKEFILE = SOURCE / "Makefile.am"

def fail(message: str, findings: list[str]) -> None:
    findings.append(message)

def main() -> int:
    findings: list[str] = []
    sources = sorted(p for p in SOURCE.rglob("*.f90")
                     if "python_diagnostics" not in p.parts)
    modules: dict[str, Path] = {}
    contents: dict[Path, str] = {}

    for path in sources:
        content = path.read_text(encoding="utf-8")
        # Strip Fortran comments without accidentally treating strings
        # as comments (the module and USE statements contain no quoted text).
        code = "\n".join(line.split("!")[0] for line in content.splitlines())
        contents[path] = code
        for m in re.finditer(
            r"(?im)^\s*module\s+(?!procedure\b|subroutine\b|function\b)(\w+)",
            code,
        ):
            name = m.group(1).lower()
            if name in modules and modules[name] != path:
                fail(f"duplicate MODULE {name}: {modules[name]} and {path}", findings)
            modules[name] = path

    make = MAKEFILE.read_text(encoding="utf-8")
    unfolded = re.sub(r"\\\s*\n", " ", make)
    match = re.search(r"(?ms)^IMEX_SfloW2D_SOURCES\s*=\s*(.*?)(?=\n\s*\n)", make)
    if match is None:
        fail("missing IMEX_SfloW2D_SOURCES in Makefile.am", findings)
        registered: set[str] = set()
    else:
        registered = set(re.findall(r"[\w/.-]+\.f90", match.group(1)))
    for path in sources:
        relative = path.relative_to(SOURCE).as_posix()
        if relative not in registered:
            fail(f"source omitted from IMEX_SfloW2D_SOURCES: {relative}", findings)

    deps: dict[str, set[str]] = defaultdict(set)
    for line in unfolded.splitlines():
        m = re.match(r"^([\w/.-]+)\.\$\(OBJEXT\)\s*:\s*(.*)", line)
        if m:
            for dep in re.finditer(r"([\w/.-]+)\.\$\(OBJEXT\)|([\w/.-]+\.inc)\b", m.group(2)):
                deps[m.group(1)].add(dep.group(1) or dep.group(2))

    for path, content in contents.items():
        rel = path.relative_to(SOURCE).as_posix()
        object_name = rel[:-4]
        for m in re.finditer(
            r"(?im)^\s*use\b\s*(?:,\s*(?:non_)?intrinsic\s*::\s*|::\s*)?(\w+)",
            content,
        ):
            used_module = m.group(1).lower()
            provider = modules.get(used_module)
            if provider is None or provider == path:
                continue  # external module or self-contained source
            expected = provider.relative_to(SOURCE).as_posix()[:-4]
            if expected not in deps[object_name]:
                fail(f"missing USE prerequisite: {object_name} -> {expected} ({used_module})", findings)

        for m in re.finditer(r"""(?im)^\s*include\s*['"]([^'"]+\.inc)['"]""", content):
            included = m.group(1)
            include_file = path.parent / included
            if not include_file.is_file():
                fail(f"missing INCLUDE file: {rel}: {included}", findings)
                continue
            expected = include_file.relative_to(SOURCE).as_posix()
            if expected not in deps[object_name]:
                fail(f"missing INCLUDE prerequisite: {object_name} -> {expected}", findings)

        category = path.relative_to(SOURCE).parts
        if category[0] == "physics":
            basename = path.stem
            required = (
                len(category) >= 3
                and category[1] in {"rheology_laws", "mass_exchange",
                                   "mass_injection_terms", "pore_pressure",
                                   "thermal_laws"}
                and basename not in {
                    "rheology_analytic_utils_2d", "rheology_parameters_2d",
                    "mass_injection_coefficients_2d"}
            ) or basename == "source_terms_2d"
            if required:
                names = set(x.lower() for x in re.findall(
                    r"(?im)^\s*(?:pure\s+)?(?:subroutine|interface)\s+(\w+)",
                    content,
                ))
                for suffix in ("explicit", "semi_implicit", "implicit", "analytic"):
                    if not any(n.endswith("_" + suffix) for n in names):
                        fail(f"source API missing {suffix}: {rel}", findings)
            if re.search(r"(?i)\b(?:q|qp)\s*\(\s*[jk]\s*,\s*[jk]\s*,", content):
                fail(f"potential direct global-state indexing in physics: {rel}", findings)

    state: dict[str, int] = {}
    stack: list[str] = []
    cycles: set[tuple[str, ...]] = set()

    def visit(node: str) -> None:
        if state.get(node) == 2:
            return
        if state.get(node) == 1:
            pos = stack.index(node)
            cycles.add(tuple(stack[pos:] + [node]))
            return
        state[node] = 1
        stack.append(node)
        for dep in sorted(deps.get(node, ())):
            if not dep.endswith(".inc"):
                visit(dep)
        stack.pop()
        state[node] = 2

    for node in sorted(deps):
        visit(node)
    for cycle in sorted(cycles):
        fail("dependency cycle: " + " -> ".join(cycle), findings)

    print(f"Fortran audit: {len(sources)} sources, {len(modules)} modules, "
          f"{sum(len(x) for x in deps.values())} declared prerequisites")
    if findings:
        print("\n".join("ERROR: " + x for x in findings), file=sys.stderr)
        return 1
    print("PASS: source registry, four-interface presence, direct USE/INCLUDE "
          "prerequisites, and dependency acyclicity")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
