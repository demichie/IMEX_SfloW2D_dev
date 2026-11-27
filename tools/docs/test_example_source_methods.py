"""Check that public examples make their source integration schedule explicit."""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "EXAMPLES"


def namelist(text, name):
    match = re.search(r"^\s*&" + re.escape(name) + r"\b(.*?)^\s*/\s*$",
                      text, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    return match.group(1) if match else None


def setting(text, name):
    return re.search(r"^\s*" + re.escape(name) + r"\s*=\s*([^,\n]+)",
                     text, re.IGNORECASE | re.MULTILINE)


class ExampleSourceMethodsTests(unittest.TestCase):
    def assert_schedule(self, section, prefix):
        self.assertIsNotNone(section, f"Missing namelist for {prefix}")
        for suffix in ("METHOD", "GROUP"):
            self.assertIsNotNone(
                setting(section, prefix + suffix),
                f"Missing explicit {prefix + suffix}",
            )

    def test_examples_declare_source_schedules(self):
        inputs = sorted(EXAMPLES.rglob("IMEX_SfloW2D.template"))
        inputs += sorted(EXAMPLES.rglob("IMEX_SfloW2D.inp"))
        self.assertTrue(inputs)
        for path in inputs:
            with self.subTest(path=path.relative_to(ROOT)):
                text = path.read_text(encoding="utf-8")
                newrun = namelist(text, "NEWRUN_PARAMETERS") or ""
                self.assert_schedule(namelist(text, "TEMPERATURE_PARAMETERS")
                                     or namelist(text, "TEMPERATURE_L1_PARAMETERS"),
                                     "THERMAL%THERMAL_" if "TEMPERATURE_PARAMETERS" in text.upper()
                                     else "THERMAL_L1%THERMAL_")
                if namelist(text, "TEMPERATURE_L2_PARAMETERS") is not None:
                    self.assert_schedule(namelist(text, "TEMPERATURE_L2_PARAMETERS"),
                                         "THERMAL_L2%THERMAL_")
                self.assert_schedule(namelist(text, "MASS_EXCHANGE_PARAMETERS")
                                     or namelist(text, "MASS_EXCHANGE_L1_PARAMETERS"),
                                     "MASS_EXCHANGE%" if namelist(text, "MASS_EXCHANGE_PARAMETERS") is not None
                                     else "MASS_EXCHANGE_L1%")
                if namelist(text, "MASS_EXCHANGE_L2_PARAMETERS") is not None:
                    self.assert_schedule(namelist(text, "MASS_EXCHANGE_L2_PARAMETERS"), "MASS_EXCHANGE_L2%")
                self.assert_schedule(namelist(text, "EXPL_TERMS_PARAMETERS"), "GEOMETRY_CURVATURE_")
                if re.search(r"^\s*RHEOLOGY_FLAG\s*=\s*T\b", newrun, re.IGNORECASE | re.MULTILINE):
                    self.assert_schedule(namelist(text, "RHEOLOGY_PARAMETERS")
                                         or namelist(text, "RHEOLOGY_L1_PARAMETERS"),
                                         "RHEOLOGY%RHEOLOGY_" if namelist(text, "RHEOLOGY_PARAMETERS") is not None
                                         else "RHEOLOGY_L1%RHEOLOGY_")
                if re.search(r"^\s*RHEOLOGY_FLAG\s*=\s*T\b", newrun,
                             re.IGNORECASE | re.MULTILINE):
                    l2 = namelist(text, "RHEOLOGY_L2_PARAMETERS")
                    if l2 is not None and not re.search(
                            r"^\s*RHEOLOGY_L2%RHEOLOGY_MODEL\s*=\s*['\"]none['\"]",
                            l2, re.IGNORECASE | re.MULTILINE):
                        self.assert_schedule(l2, "RHEOLOGY_L2%RHEOLOGY_")
                if re.search(r"^\s*RADIAL_SOURCE_FLAG\s*=\s*T\b", newrun, re.IGNORECASE | re.MULTILINE):
                    self.assert_schedule(namelist(text, "RADIAL_SOURCE_PARAMETERS"), "RADIAL_INJECTION_")
                if re.search(r"^\s*BOTTOM_RADIAL_SOURCE_FLAG\s*=\s*T\b", newrun, re.IGNORECASE | re.MULTILINE):
                    self.assert_schedule(namelist(text, "RADIAL_SOURCE_PARAMETERS"), "BOTTOM_CIRCULAR_")
                if re.search(r"^\s*BOTTOM_FISSURAL_SOURCE_FLAG\s*=\s*T\b", newrun, re.IGNORECASE | re.MULTILINE):
                    self.assert_schedule(namelist(text, "FISSURAL_SOURCE_PARAMETERS"), "BOTTOM_FISSURAL_")
                if re.search(r"^\s*PORE_PRESSURE_FLAG\s*=\s*T\b", newrun, re.IGNORECASE | re.MULTILINE):
                    self.assert_schedule(namelist(text, "PORE_PRESSURE_PARAMETERS"), "PORE_PRESSURE%PORE_PRESSURE_")

    def test_runnable_examples_document_source_integration(self):
        inputs = sorted(EXAMPLES.rglob("IMEX_SfloW2D.template"))
        inputs += sorted(EXAMPLES.rglob("IMEX_SfloW2D.inp"))
        for path in inputs:
            example_dir = path.parent
            if example_dir.name == "templatedir":
                example_dir = example_dir.parent
            readme = example_dir / "README.md"
            with self.subTest(path=readme.relative_to(ROOT)):
                self.assertTrue(readme.is_file())
                body = readme.read_text(encoding="utf-8-sig")
                self.assertIn("## Source integration", body)
                section = body.split("## Source integration", 1)[1].split("\n## ", 1)[0]
                self.assertIn("| Group | Scope | Method | Status | Role |", section)
                for heading in ("**Group**", "**Scope**", "**Method**",
                                "**Status**", "**Role**"):
                    self.assertIn(heading, section)
                group_ids = [int(value) for value in re.findall(
                    r"^\\|\\s*(\\d+)\\s*\\|", section, re.MULTILINE)]
                self.assertTrue(group_ids, "Missing source group rows")
                self.assertEqual(group_ids, sorted(group_ids))


if __name__ == "__main__":
    unittest.main()
