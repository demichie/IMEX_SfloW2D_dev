"""Prevent development-only content and broken links in the public user manual."""

import re
import unittest
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[2]
SITE = ROOT / "docs" / "website"
EXPECTED = {
    "index.md",
    "model/index.md",
    "model/numerics.md",
    "model/rheology.md",
    "model/mass-exchange.md",
    "model/thermal.md",
    "model/additional-processes.md",
    "model/applications.md",
    "model/two-layer.md",
    "model/coupling.md",
    "model/interlayer-processes.md",
    "model/hyperbolicity.md",
    "getting-started/index.md",
    "getting-started/initial-conditions.md",
    "getting-started/netcdf-fields.md",
    "getting-started/postprocessing.md",
    "references/index.md",
}
LINK = re.compile(r"!?\[[^\]\n]*\]\(([^)\n]+)\)")
DEV_TERMS = re.compile(r"\b(?:N7|TL13[ab]?(?:-[ABC])?|workflow_dispatch)\b", re.I)


class PublicWebsiteTests(unittest.TestCase):
    def test_user_documentation_is_self_contained(self):
        for rel in EXPECTED:
            with self.subTest(page=rel):
                self.assertTrue((SITE / rel).is_file())

        self.assertFalse((SITE / "development").exists())
        for path in sorted(SITE.rglob("*.md")):
            text = path.read_text(encoding="utf-8")
            with self.subTest(page=str(path.relative_to(SITE))):
                self.assertIsNone(DEV_TERMS.search(text))
                self.assertNotIn("/blob/master/docs/", text)
                self.assertNotIn("docs/acceptance", text)

    def test_development_disclaimer_is_prominent(self):
        website = (SITE / "index.md").read_text(encoding="utf-8")
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        for page in (website, readme):
            self.assertIn("under active development", page)
            self.assertIn("not intended for production use", page)
        self.assertIn('!!! danger "Development version', website)
        self.assertIn("> [!WARNING]", readme)

    def test_project_logo_sources_and_configuration(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        config = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
        builder = (ROOT / "tools" / "docs" / "build_examples.py").read_text(encoding="utf-8")
        self.assertIn("![IMEX-SfloW2D logo](docs/IMEX-SfloW2D.png)", readme)
        self.assertIn("logo: assets/images/IMEX-SfloW2D.svg", config)
        self.assertIn("favicon: assets/images/IMEX-SfloW2D.png", config)
        for filename in ("IMEX-SfloW2D.svg", "IMEX-SfloW2D.png"):
            with self.subTest(filename=filename):
                self.assertTrue((ROOT / "docs" / filename).is_file())
                self.assertIn('repo / "docs" / filename', builder)

    def test_internal_markdown_links_resolve(self):
        for path in sorted(SITE.rglob("*.md")):
            content = path.read_text(encoding="utf-8")
            for match in LINK.finditer(content):
                raw = match.group(1).strip().strip("<>")
                url = urlsplit(raw)
                if not url.path or url.scheme or raw.startswith("//"):
                    continue
                target = (path.parent / unquote(url.path)).resolve()
                with self.subTest(page=str(path.relative_to(SITE)), link=raw):
                    self.assertTrue(
                        target.is_relative_to(SITE.resolve()),
                        f"Local site link leaves docs/website: {raw}",
                    )
                    if target.is_file() or (target / "index.md").is_file():
                        continue
                    relative = target.relative_to(SITE.resolve())
                    if relative.parts and relative.parts[0] == "examples":
                        # Catalogue pages are generated from the corresponding
                        # EXAMPLES/**/README.md by build_examples.py.
                        if relative == Path("examples/index.md"):
                            continue
                        self.assertEqual(relative.name, "index.md")
                        example = ROOT / "EXAMPLES" / Path(*relative.parts[1:-1])
                        self.assertTrue((example / "README.md").is_file())
                        continue
                    self.fail(f"Unresolved website link: {raw} -> {relative}")

    def test_documented_example_solver_paths(self):
        for readme in sorted((ROOT / "EXAMPLES").glob("EXAMPLE_*/**/README.md")):
            with self.subTest(example=str(readme.relative_to(ROOT))):
                body = readme.read_text(encoding="utf-8-sig")
                self.assertNotIn("../../bin/IMEX_SfloW2D", body)
                self.assertNotIn("../../../bin/IMEX_SfloW2D", body)

    def test_references_are_bibliographic_not_developer_navigation(self):
        page = (SITE / "references" / "index.md").read_text(encoding="utf-8")
        self.assertIn("CITATION.cff", page)
        self.assertIn("de' Michieli Vitturi", page)
        self.assertNotIn("Numerical architecture", page)
        self.assertNotIn("Scientific acceptance records", page)


if __name__ == "__main__":
    unittest.main()
