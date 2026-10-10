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


def instrument(source):
    """Insert I/O immediately before the single production IMEX advance call."""
    marker = "      CALL simulation%time_integration%advance("
    if source.count(marker) != 1:
        raise AssertionError("unexpected main signature: exactly one advance hook required")
    hook = """      BLOCK
         INTEGER :: accepted_time_unit
         OPEN(NEWUNIT=accepted_time_unit,FILE='accepted_times.bin', &
              ACCESS='STREAM',FORM='UNFORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(accepted_time_unit) simulation%runtime%t, simulation%runtime%dt
         CLOSE(accepted_time_unit)
      END BLOCK
"""
    return source.replace(marker, hook+marker)


def configured_flags(directory):
    """Use the actual configured FCFLAGS; standalone strict builds have no Makefile."""
    makefile = directory/"Makefile"
    if makefile.exists():
        match = re.search(r"^FCFLAGS\s*=\s*(.*)$", makefile.read_text(), re.M)
        if match:
            return shlex.split(match[1])
    return ["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]


def build(executable, directory, repository):
    """Link an observational main against the original numerical object files."""
    executable = Path(os.environ.get("IMEX_AUDIT_SOLVER", str(executable))).resolve()
    directory.mkdir()
    source = repository/"src"
    make = (source/"Makefile.am").read_text()
    section = make.split("IMEX_SfloW2D_SOURCES =", 1)[1].split("\n\n", 1)[0]
    names = section.replace("\\", " ").split()
    numerical = [name for name in names if name != "IMEX_SfloW2D.f90"]
    objects = [executable.parent/Path(name).with_suffix(".o") for name in numerical]
    flags = configured_flags(executable.parent)+["-fopenmp"]
    includes = shlex.split(subprocess.check_output(["nf-config", "--fflags"], text=True))
    libraries = ["-llapack"]+shlex.split(subprocess.check_output(["nf-config", "--flibs"], text=True))
    libraries += shlex.split(subprocess.check_output(["nc-config", "--libs"], text=True))
    if not all(path.exists() for path in objects):
        # A supplied binary may not retain its objects. Compile the same
        # sources; field-byte equivalence below remains mandatory, not assumed.
        objects = []
        for name in numerical:
            subprocess.run(["gfortran", *flags, *includes, "-c", str(source/name)],
                           cwd=directory, check=True)
            objects.append(directory/Path(name).with_suffix(".o"))
        module_directory = directory
    else:
        module_directory = executable.parent
    main = directory/"observed_main.f90"
    main.write_text(instrument((source/"IMEX_SfloW2D.f90").read_text()))
    output = directory/"observe_accepted_times"
    subprocess.run(["gfortran", *flags, *includes, "-I"+str(module_directory),
                    str(main), *(str(path) for path in objects), *libraries,
                    "-o", str(output)], cwd=directory, check=True)
    print(output, flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("executable", type=Path)
    cli.add_argument("directory", type=Path)
    cli.add_argument("repository", type=Path)
    args = cli.parse_args()
    build(args.executable.resolve(), args.directory.resolve(), args.repository.resolve())
