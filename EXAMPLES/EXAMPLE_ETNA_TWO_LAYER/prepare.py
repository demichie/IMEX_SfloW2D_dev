#!/usr/bin/env python3
"""Prepare the Etna L1-only initial state and a strictly dry L2.

The original Etna ellipsoid script and the single-layer solver are authoritative:
no constitutive formulas or interpolation rules are reimplemented in Python.
The production N_LAYERS=2 gate is closed unless validation is explicitly opted in.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile


EXAMPLE = Path(__file__).resolve().parent
ETNA = EXAMPLE.parent / "EXAMPLE_ETNA"
ARCHIVE = ETNA / "DEM" / "Etna2014_crop.zip"
DEM_MEMBER = "Etna2014_crop.asc"
NX, NY, N_VARS = 250, 150, 5
FIELD_WIDTH = 20
SEED_NAME = "etna_initial_seed"


def replace_single_line(text, key, value):
    """Replace one whole NAMELIST assignment without editing other fields."""
    import re

    expression = re.compile(r"(?im)^([ \t]*)" + re.escape(key) + r"[ \t]*=.*$")
    result, count = expression.subn(lambda match: match.group(1) + key + "=" + value + ",", text)
    if count != 1:
        raise ValueError(f"Expected exactly one {key} in the original Etna input, found {count}")
    return result


def seed_input():
    """Request only the original single-layer initial conservative snapshot."""
    text = (ETNA / "IMEX_SfloW2D.inp").read_text(encoding="utf-8")
    for name, value in (
        ("RUN_NAME", "'" + SEED_NAME + "'"),
        ("T_END", "0.0"),
        ("OUTPUT_CONS_FLAG", "T"),
        ("OUTPUT_ESRI_FLAG", "F"),
        ("OUTPUT_NETCDF_FLAG", "F"),
        ("OUTPUT_RUNOUT_FLAG", "F"),
    ):
        text = replace_single_line(text, name, value)
    return text


def prepare_terrain(directory):
    """Extract the one canonical DEM and run the unchanged Etna preprocessor."""
    dem = directory / "DEM"
    dem.mkdir(parents=True, exist_ok=True)
    original = dem / DEM_MEMBER
    new_dem = dem / "new_dem.asc"
    initial_volume = dem / "init_volume.asc"
    if original.is_symlink() or new_dem.is_symlink() or initial_volume.is_symlink():
        raise ValueError("Refusing symlinked DEM files")
    if not original.exists():
        with zipfile.ZipFile(ARCHIVE) as archive:
            with archive.open(DEM_MEMBER) as source, original.open("xb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
    if new_dem.exists() != initial_volume.exists():
        raise ValueError("Incomplete Etna preprocessing: remove or restore both generated rasters")
    if not new_dem.exists():
        subprocess.run([sys.executable, str(ETNA / "create_input_ellipsoid.py")],
                       cwd=directory, check=True)


def write_layer_snapshots(source, l1, l2, nx=NX, ny=NY, n_vars=N_VARS):
    """Copy the solver's initial L1 snapshot exactly; make every L2 field zero.

    The fixed-width records are identical to the solver's Fortran
    (2e20.12,100(e20.12)) output, including cell order and row separators.
    Both output files are published only after complete validation.
    """
    if l1.exists() or l2.exists() or l1.is_symlink() or l2.is_symlink():
        raise FileExistsError("Initial layer snapshots already exist; refusing to overwrite")
    l1.parent.mkdir(parents=True, exist_ok=True)
    l2.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    rows = 0
    trailing_separators = 0
    solid_mass = 0.0
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="ascii", newline="\n",
                                         dir=l2.parent, prefix=".l2-", delete=False) as output:
            temp_path = Path(output.name)
            with source.open("r", encoding="ascii") as original:
                for line in original:
                    if not line.strip():
                        if count % nx != 0:
                            raise ValueError("Blank separator inside a grid row")
                        if rows < ny:
                            rows += 1
                        else:
                            # The Fortran writer appends two blank records
                            # after the final grid-row separator.
                            if count != nx * ny or trailing_separators >= 2:
                                raise ValueError("Unexpected trailing snapshot separator")
                            trailing_separators += 1
                        output.write("\n")
                        continue
                    if rows != count // nx:
                        raise ValueError("Missing grid-row separator")
                    fields = line.rstrip("\r\n")
                    if len(fields) != FIELD_WIDTH * (n_vars + 2):
                        raise ValueError("Invalid fixed-width conservative snapshot record")
                    try:
                        values = [float(fields[i:i + FIELD_WIDTH].replace("D", "E"))
                                  for i in range(0, len(fields), FIELD_WIDTH)]
                    except ValueError as exc:
                        raise ValueError("Non-numeric conservative snapshot record") from exc
                    from math import isfinite
                    if not all(map(isfinite, values)):
                        raise ValueError("Non-finite initial conservative state")
                    if values[2] < 0.0:
                        raise ValueError("Negative L1 cell mass")
                    solid_mass += values[2]
                    output.write(fields[:2 * FIELD_WIDTH])
                    output.write("".join(f"{0.0:20.12E}" for _ in range(n_vars)))
                    output.write("\n")
                    count += 1
        if count != nx * ny or rows != ny or solid_mass <= 0.0:
            raise ValueError(
                f"Invalid Etna initial snapshot: {count} cells, {rows} rows, L1 mass={solid_mass}"
            )
        shutil.copyfile(source, l1)
        os.replace(temp_path, l2)
        temp_path = None
    except Exception:
        if l1.exists():
            l1.unlink()
        raise
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def prepare(directory, solver):
    """Create the two-layer input data without invoking the disabled two-layer driver."""
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    initial = directory / "INITIAL_CONDITIONS"
    l1 = initial / "etna_L1_0000.q_2d"
    l2 = initial / "etna_L2_0000.q_2d"
    if l1.exists() or l2.exists():
        raise FileExistsError("Initial snapshots exist; remove both explicitly to regenerate")
    prepare_terrain(directory)
    binary = solver.resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"Build the single-layer solver first: {binary}")
    with tempfile.TemporaryDirectory(prefix="etna-initial-") as temp_name:
        temp = Path(temp_name)
        (temp / "DEM").symlink_to((directory / "DEM").resolve(), target_is_directory=True)
        (temp / "IMEX_SfloW2D.inp").write_text(seed_input(), encoding="utf-8")
        subprocess.run([str(binary)], cwd=temp, check=True,
                       env={**os.environ, "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"})
        matches = list((temp / "OUTPUT" / SEED_NAME).rglob(f"{SEED_NAME}_0000.q_2d"))
        if len(matches) != 1:
            raise RuntimeError(f"Expected one t=0 Etna conservative snapshot, found {len(matches)}")
        write_layer_snapshots(matches[0], l1, l2)
    print("L1: original Etna conservative state:", l1)
    print("L2: zero mass, momentum, energy and solid content:", l2)
    print("N_LAYERS=2 requires IMEX_TWO_LAYER_VALIDATION=1 for controlled testing.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=EXAMPLE,
                        help="Target example directory (default: alongside prepare.py)")
    parser.add_argument("--solver", type=Path,
                        default=EXAMPLE.parent.parent / "bin" / "IMEX_SfloW2D",
                        help="Path to the compiled single-layer solver")
    args = parser.parse_args()
    prepare(args.directory, args.solver)


if __name__ == "__main__":
    main()
