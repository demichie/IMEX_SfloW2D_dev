#!/usr/bin/env python3
"""Add a small Ocean Blue & Teal override to generated Doxygen HTML."""
import sys
from pathlib import Path

root = Path(sys.argv[1])
if not (root / "index.html").exists():
    raise SystemExit("Doxygen index.html missing")
css = """/* IMEX_SfloW2D Ocean Blue & Teal */
:root { --primary-color: #153B5B; --primary-dark-color: #102D45;
--primary-light-color: #E5F1F3; --primary-lighter-color: #E5F1F3;
--link-color: #247D98; --search-result-link-color: #247D98; }
a, a:visited { color: #247D98; }
#top, #titlearea { border-color: #247D98; }
"""
(root / "imex-theme.css").write_text(css, encoding="utf-8")
for page in root.rglob("*.html"):
    raw = page.read_text(encoding="utf-8", errors="replace")
    if "</head>" in raw:
        # Doxygen output HTML normally lives in the top-level doxygen directory.
        relative_css = "../" * (len(page.relative_to(root).parts) - 1) + "imex-theme.css"
        raw = raw.replace("</head>", f'<link rel="stylesheet" href="{relative_css}" />\n</head>', 1)
        page.write_text(raw, encoding="utf-8")
print("Doxygen styling applied")
