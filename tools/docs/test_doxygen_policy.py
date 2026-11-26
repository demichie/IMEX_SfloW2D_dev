"""Verify that Doxygen publication enforces source documentation contracts."""

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def option(name):
    text = (ROOT / "template.doxy").read_text(encoding="utf-8")
    matches = re.findall(rf"(?m)^{re.escape(name)}[ \t]*=[ \t]*(.*?)[ \t]*$", text)
    if len(matches) != 1:
        raise ValueError(f"Expected one Doxygen setting for {name}, found {matches!r}")
    return matches[0]


class DoxygenPolicyTests(unittest.TestCase):
    def test_strict_documentation_policy(self):
        expected = {
            "EXTRACT_ALL": "NO",
            "EXTRACT_PRIVATE": "YES",
            "WARNINGS": "YES",
            "WARN_IF_DOC_ERROR": "YES",
            "WARN_IF_INCOMPLETE_DOC": "YES",
            "WARN_NO_PARAMDOC": "YES",
            "WARN_AS_ERROR": "FAIL_ON_WARNINGS",
        }
        for name, value in expected.items():
            with self.subTest(setting=name):
                self.assertEqual(option(name), value)

    def test_only_one_authoritative_main_page(self):
        self.assertEqual(option("USE_MDFILE_AS_MAINPAGE"), "")
        program = (ROOT / "src" / "IMEX_SfloW2D.f90").read_text(encoding="utf-8")
        self.assertEqual(len(re.findall(r"(?m)^\s*!>\s*[@\\]mainpage\b", program)), 1)

    def test_warning_log_is_retained(self):
        self.assertEqual(option("WARN_LOGFILE"), "build/doxygen.warnings")


if __name__ == "__main__":
    unittest.main()
