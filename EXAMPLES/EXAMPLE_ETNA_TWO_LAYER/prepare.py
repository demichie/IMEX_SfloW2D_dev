#!/usr/bin/env python3
"""Generate Etna terrain and thickness rasters without running the solver.

The two-layer Fortran solver reads layer-specific primitive input directly.
The original Etna DEM must be extracted manually before preprocessing.
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


EXAMPLE = Path(__file__).resolve().parent
DEM_MEMBER = "Etna2014_crop.asc"
NX, NY, N_VARS = 250, 150, 5
FIELD_WIDTH = 20


def replace_single_line(text, key, value):
    """Replace one whole NAMELIST assignment without editing other fields."""
    import re

    expression = re.compile(r"(?im)^([ \t]*)" + re.escape(key) + r"[ \t]*=.*$")
    result, count = expression.subn(lambda match: match.group(1) + key + "=" + value + ",", text)
    if count != 1:
        raise ValueError(f"Expected exactly one {key} in the original Etna input, found {count}")
    return result





def prepare_terrain(directory):
    """Generate terrain rasters from the manually extracted original DEM."""
    dem_dir = directory / "DEM"
    original_dem = dem_dir / DEM_MEMBER
    if not original_dem.is_file():
        raise FileNotFoundError(
            f"Missing {original_dem}. Extract the original DEM manually before preprocessing."
        )

    output_files = (dem_dir / "new_dem.asc", dem_dir / "init_volume.asc",
                    dem_dir / "ellipsoid.asc")
    if all(path.is_file() for path in output_files):
        return

    subprocess.run([sys.executable, str(EXAMPLE / "create_input_ellipsoid.py")],
                   cwd=directory, check=True)


def prepare(directory):
    """Prepare the terrain rasters only; never run the Fortran executable."""
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if directory != EXAMPLE:
        original_dem = EXAMPLE / "DEM" / DEM_MEMBER
        target_dem = directory / "DEM" / DEM_MEMBER
        if not target_dem.is_file():
            if not original_dem.is_file():
                raise FileNotFoundError(
                    f"Missing {original_dem}. Extract the DEM manually first."
                )
            target_dem.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original_dem, target_dem)
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
