"""Check exact Etna lava input/DEM preparation, including non-overwrite controls."""

import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile


REPO = Path(__file__).resolve().parents[2]
EXAMPLE = REPO / "EXAMPLES/EXAMPLE_ETNA_LAVA"
SPEC = importlib.util.spec_from_file_location("etna_lava_prepare", EXAMPLE / "prepare.py")
prepare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prepare)


class PreparationTests(unittest.TestCase):
    def archive(self, directory, member=prepare.MEMBER):
        """Create small controlled ZIP data, unrelated to the production raster."""
        path = directory / "test.zip"
        content = b"original unprocessed DEM\n"
        with zipfile.ZipFile(path, "w") as output:
            output.writestr(member, content)
        return path, prepare.digest(path), hashlib.sha256(content).hexdigest()

    def test_supplied_input_is_frozen_without_parameter_changes(self):
        fixture = json.loads((EXAMPLE / "fixture.json").read_text())
        self.assertEqual(prepare.digest(EXAMPLE / "IMEX_SfloW2D.inp"), fixture["input"]["sha256"])
        self.assertFalse(fixture["input"]["parameters_changed_on_import"])
        self.assertEqual(fixture["N9"]["selected_long_examples"],
                         ["EXAMPLE_ETNA", "EXAMPLE_ETNA_LAVA"])
        self.assertFalse(fixture["N9"]["scientific_baseline_approved"])

    def test_original_archive_and_raster_cover_the_simulation_grid(self):
        fixture = json.loads((EXAMPLE / "fixture.json").read_text())
        with tempfile.TemporaryDirectory(prefix="imex-etna-lava-setup-") as temporary:
            directory = Path(temporary)
            raster = prepare.prepare(directory)
            self.assertEqual(prepare.digest(raster), fixture["DEM"]["member_sha256"])
            self.assertEqual(raster.stat().st_size, fixture["DEM"]["bytes"])
            self.assertEqual((directory / "IMEX_SfloW2D.inp").read_bytes(),
                             (EXAMPLE / "IMEX_SfloW2D.inp").read_bytes())
            with raster.open() as stream:
                header = {key.lower(): float(value) for key, value in
                          (stream.readline().split() for _ in range(6))}
            self.assertLessEqual(header["xllcorner"], 499813)
            self.assertGreaterEqual(header["xllcorner"] + header["ncols"] * header["cellsize"],
                                    499813 + 250 * 10)
            self.assertLessEqual(header["yllcorner"], 4176700)
            self.assertGreaterEqual(header["yllcorner"] + header["nrows"] * header["cellsize"],
                                    4176700 + 150 * 10)
            before = raster.stat().st_mtime_ns
            prepare.prepare(directory)
            self.assertEqual(raster.stat().st_mtime_ns, before)

    def test_small_archive_identity_and_idempotence(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            archive, archive_hash, member_hash = self.archive(directory)
            destination = directory / "DEM" / prepare.MEMBER
            prepare.install_dem(archive, destination, archive_hash, member_hash)
            before = destination.read_bytes()
            prepare.install_dem(archive, destination, archive_hash, member_hash)
            self.assertEqual(destination.read_bytes(), before)

    def test_existing_different_dem_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            archive, archive_hash, member_hash = self.archive(directory)
            destination = directory / prepare.MEMBER
            destination.write_bytes(b"user data")
            with self.assertRaisesRegex(ValueError, "will not be overwritten"):
                prepare.install_dem(archive, destination, archive_hash, member_hash)
            self.assertEqual(destination.read_bytes(), b"user data")

    def test_bad_archive_or_member_checksum_has_no_dem_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            archive, archive_hash, member_hash = self.archive(directory)
            destination = directory / "DEM" / prepare.MEMBER
            for archive_value, member_value in (("wrong", member_hash), (archive_hash, "wrong")):
                with self.assertRaises(ValueError):
                    prepare.install_dem(archive, destination, archive_value, member_value)
                self.assertFalse(destination.exists())
                self.assertFalse(list(directory.rglob(".original-dem-*")))

    def test_named_member_does_not_extract_archive_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            archive, archive_hash, member_hash = self.archive(directory, "../unexpected.asc")
            with self.assertRaises(KeyError):
                prepare.install_dem(archive, directory / "DEM" / prepare.MEMBER,
                                    archive_hash, member_hash)
            self.assertFalse((directory / "unexpected.asc").exists())

    def test_different_input_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            input_path = directory / "IMEX_SfloW2D.inp"
            input_path.write_bytes(b"user input")
            with self.assertRaisesRegex(ValueError, "will not be overwritten"):
                prepare.prepare(directory)
            self.assertEqual(input_path.read_bytes(), b"user input")

    def test_symlink_dem_destination_is_not_followed(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            archive, archive_hash, member_hash = self.archive(directory)
            original = directory / "original.asc"
            original.write_bytes(b"user original")
            destination = directory / prepare.MEMBER
            destination.symlink_to(original)
            with self.assertRaisesRegex(ValueError, "symbolic-link"):
                prepare.install_dem(archive, destination, archive_hash, member_hash)
            self.assertEqual(original.read_bytes(), b"user original")


if __name__ == "__main__":
    unittest.main()
