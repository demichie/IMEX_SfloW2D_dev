#!/usr/bin/env python3
"""Explicitly create the Etna L1/L2 conservative initial snapshots.

This is an opt-in initialization step, separate from DEM preprocessing and
from the two-layer simulation. It runs the production Fortran executable once
with the local original single-layer input and T_END=0, without timestepping.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

from prepare import EXAMPLE, prepare_terrain, seed_input, write_layer_snapshots, SEED_NAME


def initialize(directory, solver):
    directory = directory.resolve()
    initial = directory / "INITIAL_CONDITIONS"
    l1 = initial / "etna_L1_0000.q_2d"
    l2 = initial / "etna_L2_0000.q_2d"
    if l1.exists() or l2.exists():
        raise FileExistsError("Initial snapshots exist; remove both explicitly to regenerate")
    if not (directory / "DEM" / "new_dem.asc").is_file() or not (
        directory / "DEM" / "init_volume.asc"
    ).is_file():
        raise FileNotFoundError("Run python3 prepare.py first to generate the DEM rasters")
    binary = Path(solver).resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"Compiled Fortran solver not found: {binary}")
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
    print("L1 conservative initial state:", l1)
    print("L2 strictly dry initial state:", l2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=EXAMPLE)
    parser.add_argument("--solver", type=Path, default=EXAMPLE.parent.parent / "src" / "IMEX_SfloW2D")
    args = parser.parse_args()
    initialize(args.directory, args.solver)


if __name__ == "__main__":
    main()
