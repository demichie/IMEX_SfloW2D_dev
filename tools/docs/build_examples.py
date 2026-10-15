#!/usr/bin/env python3
"""Generate MkDocs example pages from the repository EXAMPLES/**/README.md.

Usage: python3 tools/docs/build_examples.py --repo . --output build/docs-generated
"""
import argparse
import json
import re
import shutil
from pathlib import Path
from urllib.parse import quote

import yaml

REQUIRED = ("title", "summary", "category", "tags")
SECTIONS = ("Overview", "Main parameters", "How to run", "Expected output")


def parse_readme(path):
    raw = path.read_text(encoding="utf-8-sig")
    match = re.match(r"\A---\s*\n(.*?)\n---\s*\n", raw, re.S)
    if not match:
        raise ValueError(f"{path}: missing YAML front matter")
    meta = yaml.safe_load(match.group(1))
    if not isinstance(meta, dict):
        raise ValueError(f"{path}: invalid metadata")
    for key in REQUIRED:
        if not meta.get(key):
            raise ValueError(f"{path}: missing {key}")
    if not isinstance(meta["tags"], list) or not all(isinstance(t, str) for t in meta["tags"]):
        raise ValueError(f"{path}: tags must be a string list")
    body = raw[match.end():]
    missing = [s for s in SECTIONS if not re.search(r"^#{2,4}\s+" + re.escape(s) + r"\s*$", body, re.I | re.M)]
    # Do not discard existing content just because a README is still being standardized.
    return meta, body, missing


def rewrite_markdown_links(body, example_dir, repo_url, revision):
    """Preserve source-relative links/images, which would otherwise break after relocation."""
    def replacer(match):
        marker, label, target = match.groups()
        # URL fragments, absolute URL, mail and reference anchors are already resolvable.
        if target.startswith(("http://", "https://", "mailto:", "#", "/", "data:")):
            return match.group(0)
        # Preserve optional angle-bracket paths and avoid changing complex link syntax.
        url = target.split()[0].strip("<>")
        if not url or url.startswith("$"):
            return match.group(0)
        from pathlib import PurePosixPath
        pos = PurePosixPath("EXAMPLES") / PurePosixPath(example_dir.as_posix()) / url
        components = []
        for part in pos.parts:
            if part == "..":
                if components:
                    components.pop()
            elif part not in (".", ""):
                components.append(part)
        clean = "/".join(components)
        if marker == "!":
            dest = f"https://raw.githubusercontent.com/{repo_url}/{revision}/{quote(clean, safe='/')}"
        else:
            dest = f"https://github.com/{repo_url}/blob/{revision}/{quote(clean, safe='/')}"
        return f"{marker}[{label}]({dest})"
    return re.sub(r"(!?)\[([^\]\n]+)\]\(([^)\n]+)\)", replacer, body)


def generate(repo, output, revision, repo_slug):
    examples_root = repo / "EXAMPLES"
    if not examples_root.is_dir():
        raise ValueError(f"Missing {examples_root}")
    output = output.resolve()
    if output == repo.resolve() or output == examples_root.resolve():
        raise ValueError("Unsafe output directory")
    if output.exists():
        shutil.rmtree(output)
    docs = output / "docs"
    shutil.copytree(repo / "docs" / "website", docs)
    pages = docs / "examples"
    pages.mkdir(parents=True, exist_ok=True)
    records, warnings = [], []
    for readme in sorted(examples_root.rglob("README.md")):
        relative_dir = readme.parent.relative_to(examples_root)
        meta, body, missing = parse_readme(readme)
        body = rewrite_markdown_links(body, relative_dir, repo_slug, revision)
        dest = pages / relative_dir / "index.md"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(body.rstrip() + "\n\n---\n\n[Source example on GitHub](https://github.com/" + repo_slug + "/tree/" + revision + "/EXAMPLES/" + relative_dir.as_posix() + ")\n", encoding="utf-8")
        slug = relative_dir.as_posix()
        records.append((slug, meta, relative_dir))
        if missing:
            warnings.append(f"{slug}: missing sections: {', '.join(missing)}")
    if not records:
        raise ValueError("No EXAMPLES/**/README.md files discovered")
    categories = {}
    for slug, meta, rel in records:
        categories.setdefault(str(meta["category"]), []).append((slug, meta, rel))
    lines = ["# Examples", "", f"Automatically generated from **{len(records)}** example READMEs at revision `{revision[:12]}`.", "", "Every page includes a link to its original source folder. The full Markdown content, including bibliographic references and technical notes, is preserved.", ""]
    for category in sorted(categories, key=str.casefold):
        lines.extend([f"## {category}", ""])
        for slug, meta, _ in sorted(categories[category], key=lambda entry: entry[1]["title"].casefold()):
            link = f"{slug}/index.md"
            lines.extend([f"### [{meta['title']}]({link})", "", str(meta["summary"]), "", "Tags: " + ", ".join(f"`{t}`" for t in meta["tags"]), ""])
    (pages / "index.md").write_text("\n".join(lines), encoding="utf-8")
    nav = ["  - Home: index.md", "  - Getting started:", "      - Overview: getting-started/index.md", "      - Installation: getting-started/installation.md", "      - Output files: getting-started/output-files.md", "  - Model: model/index.md", "  - Examples:", "      - Catalogue: examples/index.md"]
    for category in sorted(categories, key=str.casefold):
        nav.append(f"      - {json.dumps(category, ensure_ascii=False)}:")
        for slug, meta, _ in sorted(categories[category], key=lambda entry: entry[1]["title"].casefold()):
            label = json.dumps(meta['title'], ensure_ascii=False)
            nav.append(f"          - {label}: examples/{slug}/index.md")
    nav += ["  - Development: development/index.md", "  - Doxygen: https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/", "  - References: references/index.md"]
    config = (repo / "mkdocs.yml").read_text(encoding="utf-8")
    config = config.replace("# GENERATED_NAVIGATION", "\n".join(nav))
    (output / "mkdocs.generated.yml").write_text(config, encoding="utf-8")
    (output / "example-audit.txt").write_text(f"Examples: {len(records)}\nRevision: {revision}\n" + "\n".join(warnings) + "\n", encoding="utf-8")
    print(f"Generated {len(records)} examples; {len(warnings)} content warnings")
    for warning in warnings:
        print("WARNING:", warning)
    return len(records)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--repo", type=Path, default=Path("."))
    p.add_argument("--output", type=Path, default=Path("build/docs-generated"))
    p.add_argument("--revision", default="master")
    p.add_argument("--repo-slug", default="demichie/IMEX_SfloW2D_dev")
    a = p.parse_args()
    generate(a.repo.resolve(), a.output, a.revision, a.repo_slug)
