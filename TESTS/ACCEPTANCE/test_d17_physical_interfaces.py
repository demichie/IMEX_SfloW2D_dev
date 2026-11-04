#!/usr/bin/env python3
"""Static D19 audit of production source interfaces and Lie integration contracts.

This acceptance gate intentionally inspects the current source catalog and Lie
executor; it must not depend on the removed STAGE/STEP/RK architecture.
"""
from __future__ import annotations

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
PHYSICS = SRC / "physics"
METHODS = ("explicit", "semi_implicit", "implicit", "analytic")


def compact(path: Path) -> str:
    code = path.read_text(encoding="utf-8")
    code = re.sub(r"(?m)!.*$", "", code)
    return re.sub(r"&\s*\n\s*&?", " ", code)


def has_procedure(code: str, name: str) -> bool:
    return bool(re.search(
        rf"(?im)^\s*(?:(?:pure|recursive|elemental)\s+)*"
        rf"(?:subroutine|function)\s+{re.escape(name)}\s*\(", code))


def require(ok: bool, message: str, failures: list[str]) -> None:
    if not ok:
        failures.append(message)


def main() -> int:
    failures: list[str] = []
    catalog = compact(SRC / "numerics/source_catalog_2d.f90")
    registry = compact(SRC / "numerics/source_registry_2d.f90")
    executor = compact(SRC / "numerics/lie_source_executor_2d.f90")
    assembly = compact(PHYSICS / "source_terms_2d.f90")
    integrator = compact(SRC / "numerics/time_integration_2d.f90")
    simulation = compact(SRC / "core/simulation_2d.f90")
    state = compact(SRC / "state/state_2d.f90")
    files = {path: compact(path) for path in PHYSICS.rglob("*.f90")}
    entries = re.findall(
        r"source_term_type\s*\(\s*'([^']+)'\s*,\s*"
        r"\.(true|false)\.\s*,\s*\.(true|false)\.\s*,\s*"
        r"\.(true|false)\.\s*,\s*default_method\s*=\s*"
        r"SOURCE_METHOD_(EXPLICIT|NEWTON|ANALYTIC)", catalog, re.I)
    declared_count = re.search(r"(?i)SOURCE_TERM_COUNT\s*=\s*(\d+)", catalog)
    require(declared_count is not None, "missing SOURCE_TERM_COUNT declaration", failures)
    if declared_count is not None:
        expected_count = int(declared_count.group(1))
        require(len(entries) == expected_count,
                f"expected {expected_count} source catalog entries, found {len(entries)}", failures)
    names = [e[0] for e in entries]
    require(len(set(names)) == len(names), "duplicate source identities", failures)
    expected_families = {
        "rheology": PHYSICS / "rheology_laws",
        "mass_exchange": PHYSICS / "mass_exchange",
    }
    for name, exp, newton, analytic, default in entries:
        if name.startswith("rheology_"):
            owner = expected_families["rheology"] / (name + "_2d.f90")
            base = f"eval_{name}"
            for method in METHODS:
                require(owner in files and has_procedure(files.get(owner, ""), f"{base}_{method}"),
                        f"{name}: missing {method} physical procedure", failures)
        elif name.startswith("mass_exchange_"):
            owner = expected_families["mass_exchange"] / (name + "_2d.f90")
            for method in METHODS:
                proc = f"{'solve' if method == 'analytic' else 'eval'}_{name}_{method}"
                require(owner in files and has_procedure(files.get(owner, ""), proc),
                        f"{name}: missing {method} physical procedure", failures)
        elif name == "thermal_cooling":
            owner = PHYSICS / "thermal_laws/thermal_cooling_2d.f90"
            for method in METHODS:
                require(has_procedure(files[owner], f"eval_thermal_cooling_{method}"),
                        f"thermal cooling: missing {method} procedure", failures)
        elif name == "pore_pressure_dissipation":
            owner = PHYSICS / "pore_pressure/pore_pressure_darcy_2d.f90"
            for method in METHODS:
                require(has_procedure(files[owner], f"eval_pore_pressure_darcy_{method}"),
                        f"pore pressure: missing {method} procedure", failures)
        else:
            # Injection names use a coordinator and dedicated physical owners.
            require(name in assembly or name in compact(PHYSICS / "source_configuration_2d.f90"),
                    f"unrouted source {name}", failures)
        allowed = {"EXPLICIT": exp, "NEWTON": newton, "ANALYTIC": analytic}
        require(allowed[default].lower() == "true",
                f"{name}: default METHOD {default} is not supported", failures)

    require(bool(re.search(r"interface\s+eval_source_explicit\b", assembly, re.I)),
            "missing explicit generic", failures)
    for proc in ("eval_geometry_cell_source_explicit", "eval_mass_exchange_source_explicit",
                 "assemble_lie_explicit_transaction", "eval_source_semi_implicit",
                 "eval_source_implicit", "eval_source_analytic"):
        require(has_procedure(assembly, proc), f"missing source assembler {proc}", failures)
    for fragment in ("groups(1)%id/=0", "groups(1)%method/=SOURCE_METHOD_EXPLICIT",
                     "SOURCE_METHOD_NEWTON", "source_newton_affected_equations"):
        require(fragment.lower() in registry.lower(), f"registry contract: {fragment}", failures)
    # D20 cell-major Lie execution: each worklist cell owns local candidates
    # throughout all source groups. Do not require retired full-domain group
    # snapshots (q_group/qp_group): they would duplicate storage and defeat the
    # cell-local parallel execution contract.
    for fragment in ("q_local=q(:,j,k)", "qp_local=qp(:,j,k)",
                     "qpentry=qp_original(:,j,k)", "qpentry=qp_local",
                     "q_local=q_candidate", "qp_local=qp_candidate",
                     "q(:,j,k)=q_local", "qp(:,j,k)=qp_local",
                     "CALL state%refresh_substrate_totals()"):
        require(fragment.lower() in executor.lower(), f"Lie executor contract: {fragment}", failures)
    executor_lower = executor.lower()
    ordered_fragments = ("group_loop: do group_index=1,size(executor%groups)",
                         "q(:,j,k)=q_local", "qp(:,j,k)=qp_local")
    positions = [executor_lower.find(fragment) for fragment in ordered_fragments]
    require(all(position >= 0 for position in positions) and positions == sorted(positions),
            "Lie executor contract: ordered cell-local groups before state publication", failures)
    for fragment in ("snapshot_substrate", "restore_substrate", "begin_geometry_transaction",
                     "rollback_geometry_transaction", "conservative_mass_admissible"):
        require(fragment.lower() in integrator.lower(), f"Lie rollback contract: {fragment}", failures)
    require("this%runtime%t=step_start+step_length" in simulation,
            "simulation must advance clock once", failures)
    require("CALL registry%finalize()" in compact(PHYSICS / "source_configuration_2d.f90"),
            "failed source reconfiguration may retain the previous registry", failures)
    require("SUM(q_cell(idx_add_gas_first:idx_add_gas_last))" in state,
            "missing additional-gas fraction admissibility", failures)
    require("stochastic_flag, stochastic_value, rheology_rate" in assembly,
            "explicit rheology must receive cell stochastic value", failures)

    removed = (
        "step_source_executor_2d", "eval_stage_source_explicit",
        "eval_step_source_explicit", "eval_stage_expl_terms", "eval_step_expl_terms",
    )
    for filename in SRC.rglob("*.f90"):
        production = compact(filename)
        for symbol in removed:
            if re.search(rf"\b{re.escape(symbol)}\b", production, re.I):
                failures.append(f"legacy symbol {symbol}: {filename.relative_to(ROOT)}")
    for failure in failures:
        print("FAIL:", failure, file=sys.stderr)
    print(f"D19 static source interface audit: {len(entries)} catalog entries; "
          f"{len(files)} physical files; {len(failures)} failures")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
