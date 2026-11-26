#!/usr/bin/env python3
"""Check that repository Markdown uses math syntax supported by GitHub and MkDocs."""

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FENCE = re.compile(r"^\s*(\x60{3,}|~{3,})")
CODE_SPAN = re.compile(r"(\x60+).*?\1")
UNSUPPORTED_MATH = re.compile(r"\\(?:\[|\]|\(|\))")
SKIPPED_DIRECTORIES = {".git", ".venv", "build", "OUTPUT", "__pycache__", ".tox"}


def unsupported_delimiters(content):
    """Report legacy math delimiters, ignoring fenced and inline code."""
    fence = None
    for number, line in enumerate(content.splitlines(), start=1):
        marker = FENCE.match(line)
        if marker:
            ticks = marker.group(1)
            if fence is None:
                fence = (ticks[0], len(ticks))
            elif ticks[0] == fence[0] and len(ticks) >= fence[1]:
                fence = None
            continue
        if fence is not None:
            continue
        plain_text = CODE_SPAN.sub("", line)
        if UNSUPPORTED_MATH.search(plain_text):
            yield number


class MarkdownMathTests(unittest.TestCase):
    def test_github_compatible_math_delimiters(self):
        problems = []
        for path in sorted(ROOT.rglob("*.md")):
            relative = path.relative_to(ROOT)
            if any(part in SKIPPED_DIRECTORIES for part in relative.parts):
                continue
            lines = list(unsupported_delimiters(path.read_text(encoding="utf-8-sig")))
            if lines:
                problems.append(f"{relative}: {', '.join(map(str, lines))}")
        self.assertFalse(
            problems,
            "Replace \\(...\\) with $...$ and \\[...\\] with $$...$$ outside code:\n"
            + "\n".join(problems),
        )

    def test_mkdocs_math_runtime_configured(self):
        config = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
        self.assertIn("pymdownx.arithmatex:", config)
        self.assertIn("generic: true", config)
        self.assertIn("assets/javascripts/mathjax.js", config)
        self.assertIn("mathjax@3/es5/tex-mml-chtml.js", config)
        self.assertTrue((ROOT / "docs/website/assets/javascripts/mathjax.js").is_file())


if __name__ == "__main__":
    unittest.main()
