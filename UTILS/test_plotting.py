"""Small synthetic output-contract checks for the optional plotting tools."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")


def load(name):
    location = Path(__file__).with_name(name + ".py")
    spec = importlib.util.spec_from_file_location(name, location)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def asc(path, values, nodata=-9999):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ncols 2\nnrows 2\nxllcorner 100\nyllcorner 200\n"
        "cellsize 10\nNODATA_value -9999\n"
        + "\n".join(" ".join(map(str, row)) for row in values) + "\n",
        encoding="utf-8")


class PlotTests(unittest.TestCase):
    def test_overlay_geography_and_manifest_time(self):
        overlay = load("plot_overlay")
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            asc(root / "terrain" / "dem_esri.asc", [[1, 2], [3, 4]])
            metadata = root / "metadata"
            metadata.mkdir()
            (metadata / "output_manifest.csv").write_text(
                "snapshot_index,time_seconds\n0,0.0\n1,0.7\n",
                encoding="utf-8")
            raster = root / "snapshots" / "0001_t_7E-1" / "run_0001.asc"
            asc(raster, [[2, -9999], [4, 5]])
            array, extent = overlay.read_asc(raster)
            self.assertEqual(array.shape, (2, 2))
            self.assertEqual(extent, (100, 120, 200, 220))
            self.assertTrue(overlay.manifest_times(root)[1] == 0.7)
            self.assertEqual(overlay.run(root, "snapshots"), [raster.with_suffix(".png")])
            self.assertTrue(raster.with_suffix(".png").exists())

    def test_runout_single_row_and_zero_maximum(self):
        module = load("plot_runout")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "run_runout.csv"
            path.write_text("Time (s),Runout (m),Area (m2),X,Y\n0,0,0,0,0\n",
                            encoding="utf-8")
            output, when, maximum, milestones = module.analyze(path)
            self.assertEqual((when, maximum, milestones), (0.0, 0.0, {}))
            self.assertTrue(output.exists())

    def test_probe_single_row_and_rate_rejection(self):
        module = load("plot_probes")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "run_prb_01.csv"
            path.write_text("100,200\ntime,h\n(s),(m)\n0,2\n", encoding="utf-8")
            outputs = module.plot(path)
            self.assertEqual(len(outputs), 1)
            self.assertTrue(outputs[0].exists())
            with self.assertRaisesRegex(ValueError, "Two or more"):
                module.plot(path, rate=True)


if __name__ == "__main__":
    unittest.main()
