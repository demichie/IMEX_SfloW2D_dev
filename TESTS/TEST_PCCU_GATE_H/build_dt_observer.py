"""Build a test-only main that records accepted times without numerical edits.

Reuse the tested executable's numerical objects and compiler flags when
available. Final canonical fields MUST subsequently match the uninstrumented
executable byte for byte. No source in src/ is changed by this observer.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess


def instrument(source, include_spatial=False):
    """Record accepted times before the unique production advance call.

    The CFL proposal is already saved in cfl_dt immediately after
    compute_timestep; do not duplicate or modify the numerical evaluation.
    """
    marker = re.compile(r"(?m)^(?P<indent> *)CALL simulation%advance_timestep\(")
    matches = list(marker.finditer(source))
    if len(matches) != 1:
        raise AssertionError("unexpected main signature: exactly one advance hook required")

    indent = matches[0].group("indent")
    record = "simulation%runtime%t, simulation%runtime%dt"
    if include_spatial:
        proposal = "        cfl_dt=simulation%runtime%dt"
        if source.count(proposal) != 1:
            raise AssertionError("unexpected spatial timestep proposal hook")
        record += ", cfl_dt"

    hook = (
        f"{indent}BLOCK\n"
        f"{indent}    INTEGER :: accepted_time_unit\n"
        f"{indent}    OPEN(NEWUNIT=accepted_time_unit,FILE='accepted_times.bin', &\n"
        f"{indent}        ACCESS='STREAM',FORM='UNFORMATTED',STATUS='UNKNOWN',POSITION='APPEND')\n"
        f"{indent}    WRITE(accepted_time_unit) {record}\n"
        f"{indent}    CLOSE(accepted_time_unit)\n"
        f"{indent}END BLOCK\n"
    )
    return source[:matches[0].start()] + hook + source[matches[0].start():]


def configured_flags(directory):
    """Use configured FCFLAGS, falling back to a conservative debug profile."""
    makefile = directory / "Makefile"
    if makefile.exists():
        match = re.search(r"^FCFLAGS\s*=\s*(.*)$", makefile.read_text(), re.M)
        if match:
            return shlex.split(match[1])
    compiler = shlex.split(os.environ.get("FC", "gfortran"))[0]
    if Path(compiler).name.startswith("ifx"):
        return ["-O0", "-g", "-traceback", "-check", "all"]
    return ["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]


def configured_variable(directory, name):
    """Read a configured Makefile variable when available."""
    makefile = directory / "Makefile"
    if makefile.exists():
        match = re.search(rf"^{re.escape(name)}\s*=\s*(.*)$", makefile.read_text(), re.M)
        if match:
            return shlex.split(match[1])
    return []


def build(executable, directory, repository, include_spatial=False):
    """Link an observational main against the original numerical object files."""
    executable = Path(os.environ.get("IMEX_AUDIT_SOLVER", str(executable))).resolve()
    directory.mkdir()
    source = repository / "src"
    make = (source / "Makefile.am").read_text()
    section = make.split("IMEX_SfloW2D_SOURCES =", 1)[1].split("\n\n", 1)[0]
    names = section.replace("\\", " ").split()
    numerical = [name for name in names if name != "IMEX_SfloW2D.f90"]
    objects = [executable.parent / Path(name).with_suffix(".o") for name in numerical]

    compiler = shlex.split(os.environ.get("FC", "gfortran"))
    if not compiler:
        raise ValueError("FC cannot be empty")
    openmp_flags = configured_variable(executable.parent, "OPENMP_FCFLAGS")
    if not openmp_flags:
        openmp_flags = shlex.split(os.environ.get("OPENMP_FCFLAGS", "-fopenmp"))
    flags = configured_flags(executable.parent) + openmp_flags
    includes = shlex.split(subprocess.check_output(["nf-config", "--fflags"], text=True))
    libraries = ["-llapack"] + shlex.split(subprocess.check_output(["nf-config", "--flibs"], text=True))
    libraries += shlex.split(subprocess.check_output(["nc-config", "--libs"], text=True))

    if not all(path.exists() for path in objects):
        # Never link objects from a different compiler. Rebuild when necessary.
        objects = []
        c_compiler = shlex.split(os.environ.get("CC", "cc"))
        if not c_compiler:
            raise ValueError("CC cannot be empty")
        for name in numerical:
            object_path = directory / (Path(name).stem + ".o")
            if name.endswith(".c"):
                subprocess.run([*c_compiler, "-c", str(source / name),
                                "-o", str(object_path)],
                               cwd=directory, check=True)
            else:
                subprocess.run([*compiler, *flags, *includes, "-I" + str(source / "state"),
                                "-c", str(source / name)],
                               cwd=directory, check=True)
            objects.append(object_path)
        module_directory = directory
    else:
        module_directory = executable.parent

    main = directory / "observed_main.f90"
    main.write_text(instrument((source / "IMEX_SfloW2D.f90").read_text(), include_spatial))
    output = directory / "observe_accepted_times"
    subprocess.run([*compiler, *flags, *includes, "-I" + str(module_directory),
                    str(main), *(str(path) for path in objects), *libraries,
                    "-o", str(output)], cwd=directory, check=True)
    print(output, flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("executable", type=Path)
    cli.add_argument("directory", type=Path)
    cli.add_argument("repository", type=Path)
    cli.add_argument("--include-spatial", action="store_true",
                     help="also record spatial proposals before event/growth clipping")
    args = cli.parse_args()
    build(args.executable.resolve(), args.directory.resolve(), args.repository.resolve(), args.include_spatial)
