"""Guard source distribution inventories against fixture relocation."""

import re
import unittest
from pathlib import Path

from build_archives import MAKE_CHECK_FILES


ROOT = Path(__file__).resolve().parents[2]
DIST_RE = re.compile(r"(?:^|\s)((?:TESTS|EXAMPLES)/[^\s\\]+)", re.M)


class ReleaseInventoryTests(unittest.TestCase):
    def test_required_files_exist(self):
        for path in MAKE_CHECK_FILES:
            with self.subTest(path=path):
                self.assertTrue((ROOT / path).is_file(), f"Missing release input: {path}")

    def test_required_files_are_distributed_by_automake(self):
        makefile = (ROOT / "Makefile.am").read_text(encoding="utf-8")
        dist_paths = set(DIST_RE.findall(makefile))
        for path in MAKE_CHECK_FILES:
            with self.subTest(path=path):
                self.assertIn(path, dist_paths, f"Not listed in Makefile.am EXTRA_DIST: {path}")

    def test_no_duplicate_required_entries(self):
        self.assertEqual(len(MAKE_CHECK_FILES), len(set(MAKE_CHECK_FILES)))


if __name__ == "__main__":
    unittest.main()
