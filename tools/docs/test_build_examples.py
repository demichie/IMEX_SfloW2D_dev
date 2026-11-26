"""Regression tests for MkDocs example discovery and metadata preservation."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from build_examples import generate, parse_readme


class ExampleDocumentationTests(unittest.TestCase):
    def test_mixed_readmes_and_supporting_directories(self):
        with TemporaryDirectory() as tmp:
            repo = Path(tmp)
            website = repo / "docs" / "website"
            website.mkdir(parents=True)
            (website / "index.md").write_text("# Home\n", encoding="utf-8")
            (repo / "mkdocs.yml").write_text(
                "site_name: IMEX-SfloW2D\nnav:\n# GENERATED_NAVIGATION\n",
                encoding="utf-8",
            )
            etna = repo / "EXAMPLES" / "EXAMPLE_ETNA_TWO_LAYER"
            etna.mkdir(parents=True)
            (etna / "IMEX_SfloW2D.inp").write_text("", encoding="utf-8")
            (etna / "README.md").write_text(
                "# Etna two-layer collapse\n\n"
                "A lower-layer flow and an initially dry upper current.\n\n"
                "See [the theory](../../docs/TWO_LAYER_HYPERBOLICITY.md).\n",
                encoding="utf-8",
            )
            dem = etna / "DEM"
            dem.mkdir()
            (dem / "README.md").write_text("# DEM instructions\n", encoding="utf-8")
            template = etna / "templatedir"
            template.mkdir()
            (template / "IMEX_SfloW2D.inp").write_text("", encoding="utf-8")
            (template / "README.md").write_text("# Supporting template\n", encoding="utf-8")

            bump = repo / "EXAMPLES" / "EXAMPLE_BUMP" / "SUBCRITICAL"
            bump.mkdir(parents=True)
            (bump / "IMEX_SfloW2D.template").write_text("", encoding="utf-8")
            (bump / "README.md").write_text(
                "---\ntitle: Subcritical bump\nsummary: A benchmark\n"
                "category: Analytical benchmarks\ntags: [benchmark]\n---\n"
                "# Subcritical bump\n\nBenchmark details.\n",
                encoding="utf-8",
            )
            output = repo / "build" / "docs-generated"
            self.assertEqual(generate(repo, output, "abc123", "owner/repo"), 2)
            self.assertTrue((output / "docs" / "examples" / "EXAMPLE_ETNA_TWO_LAYER" / "index.md").is_file())
            self.assertFalse((output / "docs" / "examples" / "EXAMPLE_ETNA_TWO_LAYER" / "DEM").exists())
            self.assertFalse((output / "docs" / "examples" / "EXAMPLE_ETNA_TWO_LAYER" / "templatedir").exists())
            catalogue = (output / "docs" / "examples" / "index.md").read_text(encoding="utf-8")
            self.assertIn("Etna two-layer collapse", catalogue)
            self.assertIn("Subcritical bump", catalogue)
            self.assertIn("Two-layer flows", catalogue)
            rendered = (output / "docs" / "examples" / "EXAMPLE_ETNA_TWO_LAYER" / "index.md").read_text()
            self.assertIn("https://github.com/owner/repo/blob/abc123/docs/TWO_LAYER_HYPERBOLICITY.md", rendered)
            nav = (output / "mkdocs.generated.yml").read_text()
            self.assertIn("getting-started/docker.md", nav)
            self.assertIn("model/two-layer.md", nav)
            self.assertIn("development/publishing.md", nav)

    def test_invalid_explicit_front_matter_is_rejected(self):
        with TemporaryDirectory() as tmp:
            readme = Path(tmp) / "README.md"
            readme.write_text("---\ntitle: No summary\n---\n# A\n")
            with self.assertRaisesRegex(ValueError, "missing summary"):
                parse_readme(readme)


if __name__ == "__main__":
    unittest.main()
