#!/usr/bin/env python3
"""Solver-independent tests of the standalone Etna DEM preprocessing."""
from pathlib import Path
import tempfile
import unittest

import prepare

HERE = Path(__file__).resolve().parent


class EtnaTwoLayerPreparationTest(unittest.TestCase):

    def test_terrain_requires_extracted_dem(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, "Extract the original DEM manually"):
                prepare.prepare_terrain(Path(directory))

    def test_prepare_does_not_invoke_fortran_or_create_snapshots(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dem = root / "DEM"
            dem.mkdir()
            for name in ("Etna2014_crop.asc", "new_dem.asc", "init_volume.asc", "ellipsoid.asc"):
                (dem / name).write_text("existing", encoding="ascii")
            prepare.prepare(root)
            self.assertFalse((root / "INITIAL_CONDITIONS").exists())
            self.assertFalse((root / "OUTPUT").exists())

    def test_primitive_input_is_layer_specific(self):
        text = (HERE / "IMEX_SfloW2D.inp").read_text(encoding="utf-8")
        self.assertIn("&RESTART_L1_PARAMETERS", text)
        self.assertIn("N_RESTART_FILES=1,", text)
        self.assertIn("RESTART_FILES='./DEM/init_volume.asc',", text)
        self.assertIn("&RESTART_L2_PARAMETERS", text)
        self.assertIn("N_RESTART_FILES_L2=1,", text)
        self.assertIn("RESTART_FILES_L2='NONE',", text)
        self.assertIn("&INITIAL_L1_PARAMETERS", text)
        self.assertIn("&INITIAL_L2_PARAMETERS", text)
        self.assertIn("T_INIT=373.D0", text)
        self.assertIn("SED_VOL_PERC=50.0", text)
        self.assertIn("SED_VOL_PERC_L2=0.0", text)
        self.assertNotIn(".q_2d", text)

    def test_preprocessor_files_are_example_local(self):
        self.assertTrue((HERE / "create_input_ellipsoid.py").is_file())
        self.assertTrue((HERE / "input_ellipsoid.py").is_file())


if __name__ == "__main__":
    unittest.main()
