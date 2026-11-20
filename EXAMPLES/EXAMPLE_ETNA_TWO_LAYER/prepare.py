#!/usr/bin/env python3
"""Prepare only the Etna DEM and ellipsoidal excavation rasters.

This script never launches IMEX_SfloW2D. Conservative initial snapshots
are generated separately by the explicit initialize_layers.py command.
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


EXAMPLE = Path(__file__).resolve().parent
ARCHIVE_NAME = "Etna2014_crop.zip"
SEED_INPUT = EXAMPLE / "single_layer_seed.inp"
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
    text = SEED_INPUT.read_text(encoding="utf-8")
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
    """Use the example-local ZIP and preprocessor to create the Etna rasters."""
    dem = directory / "DEM"
    dem.mkdir(parents=True, exist_ok=True)
    original = dem / DEM_MEMBER
    new_dem = dem / "new_dem.asc"
    initial_volume = dem / "init_volume.asc"
    if original.is_symlink() or new_dem.is_symlink() or initial_volume.is_symlink():
        raise ValueError("Refusing symlinked DEM files")
    if not original.exists():
        archive_path = dem / ARCHIVE_NAME
        if not archive_path.is_file():
            raise FileNotFoundError(
                f"Missing {archive_path}. Copy only the DEM ZIP from "
                "EXAMPLE_ETNA/DEM into this example's DEM directory."
            )
        with zipfile.ZipFile(archive_path) as archive:
            with archive.open(DEM_MEMBER) as source, original.open("xb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
    if new_dem.exists() != initial_volume.exists():
        raise ValueError("Incomplete Etna preprocessing: remove or restore both generated rasters")
    if not new_dem.exists():
        subprocess.run([sys.executable, str(EXAMPLE / "create_input_ellipsoid.py")],
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


def prepare(directory):
    """Prepare the DEM rasters only. Never run the Fortran executable."""
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if directory != EXAMPLE:
        local_archive = EXAMPLE / "DEM" / ARCHIVE_NAME
        target_archive = directory / "DEM" / ARCHIVE_NAME
        if not target_archive.exists():
            if not local_archive.is_file():
                raise FileNotFoundError(
                    f"Missing {local_archive}. Copy the DEM ZIP into the example DEM directory first."
                )
            target_archive.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(local_archive, target_archive)
    prepare_terrain(directory)
    print("DEM preprocessing complete. No solver has been executed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=EXAMPLE,
                        help="Target example directory (default: alongside prepare.py)")
    args = parser.parse_args()
    prepare(args.directory)

if __name__ == "__main__":
    main()
