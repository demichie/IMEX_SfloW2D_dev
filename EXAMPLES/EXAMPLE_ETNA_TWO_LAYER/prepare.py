#!/usr/bin/env python3
"""Generate the Etna DEM and initial thickness without running the solver.

The two-layer Fortran solver reads layer-specific primitive input directly.
Only the original Etna DEM ZIP must be copied from the single-layer example.
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


EXAMPLE = Path(__file__).resolve().parent
ARCHIVE_NAME = "Etna2014_crop.zip"
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
    """Extract the local DEM and generate the three terrain rasters."""
    dem_dir = directory / "DEM"
    original_dem = dem_dir / DEM_MEMBER
    output_files = (dem_dir / "new_dem.asc", dem_dir / "init_volume.asc",
                    dem_dir / "ellipsoid.asc")

    if original_dem.is_file() and all(path.is_file() for path in output_files):
        return

    if not original_dem.is_file():
        archive = dem_dir / ARCHIVE_NAME
        if not archive.is_file():
            raise FileNotFoundError(
                f"Copy only the DEM ZIP ({ARCHIVE_NAME}) into {dem_dir} before preprocessing."
            )
        dem_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive) as source:
            with source.open(DEM_MEMBER) as data, original_dem.open("wb") as target:
                shutil.copyfileobj(data, target)

    subprocess.run([sys.executable, str(EXAMPLE / "create_input_ellipsoid.py")],
                   cwd=directory, check=True)

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
