#!/usr/bin/env python3
"""Small, solver-independent regressions for Etna two-layer preparation."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("etna_two_layer_prepare", HERE / "prepare.py")
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


def record(x, y, state):
    return "".join(f"{value:20.12E}" for value in (x, y, *state)) + "\n"


class EtnaTwoLayerPreparationTest(unittest.TestCase):

    def test_l1_byte_identical_l2_exactly_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = root / "seed_0000.q_2d"
            with original.open("w", encoding="ascii") as stream:
                for j in range(2):
                    for i in range(3):
                        stream.write(record(float(i), float(j),
                                            [float(i + j + 1), 2., 3., 4., 5.]))
                    stream.write("\n")
            # Match the production Fortran writer: one separator for
            # each grid row plus two extra blank records at EOF.
            with original.open("a", encoding="ascii") as stream:
                stream.write("\n\n")
            l1 = root / "L1_0000.q_2d"
            l2 = root / "L2_0000.q_2d"
            prepare.write_layer_snapshots(original, l1, l2, nx=3, ny=2, n_vars=5)
            self.assertEqual(l1.read_bytes(), original.read_bytes())
            originals = [line for line in original.read_text().splitlines() if line.strip()]
            upper = [line for line in l2.read_text().splitlines() if line.strip()]
            self.assertEqual(len(upper), 6)
            for before, after in zip(originals, upper):
                self.assertEqual(before[:40], after[:40])
                self.assertEqual(len(after), 140)
                self.assertEqual([float(after[i:i + 20]) for i in range(40, 140, prepare.FIELD_WIDTH)],
                                 [0.0] * 5)
            with self.assertRaises(FileExistsError):
                prepare.write_layer_snapshots(original, l1, l2, nx=3, ny=2, n_vars=5)

    def test_rejects_invalid_seed_without_publishing_layers(self):
        for contents in (
            record(0., 0., [0., 0., 0., 0., 0.]) + "\n",
            record(0., 0., [1., 0., 0., 0., 0.]),
            record(0., 0., [-1., 0., 0., 0., 0.]) + "\n",
            "broken record\n\n",
            record(0., 0., [1., 0., 0., 0., 0.]) + "\n\n\n\n",
        ):
            with self.subTest(contents=contents):
                with tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    seed = root / "seed.q_2d"
                    seed.write_text(contents, encoding="ascii")
                    l1 = root / "L1.q_2d"
                    l2 = root / "L2.q_2d"
                    with self.assertRaises(ValueError):
                        prepare.write_layer_snapshots(seed, l1, l2, nx=1, ny=1, n_vars=5)
                    self.assertFalse(l1.exists())
                    self.assertFalse(l2.exists())

    def test_terrain_requires_local_dem_zip(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, "Copy only the DEM ZIP"):
                prepare.prepare_terrain(Path(directory))

    def test_seed_input_preserves_original_kinematics(self):
        text = prepare.seed_input()
        self.assertIn("T_END=0.0,", text)
        self.assertIn("OUTPUT_CONS_FLAG=T,", text)
        self.assertIn("RHEOLOGY%MU = 0.3", text)
        self.assertIn("RHEOLOGY%XI = 300.0", text)
        self.assertIn("SED_VOL_PERC = 50.0", text)
        self.assertNotIn("N_LAYERS = 2", text)
        self.assertTrue((HERE / "create_input_ellipsoid.py").is_file())
        self.assertTrue((HERE / "input_ellipsoid.py").is_file())
        self.assertTrue((HERE / "single_layer_seed.inp").is_file())


if __name__ == "__main__":
    unittest.main()
