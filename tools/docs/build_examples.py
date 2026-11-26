#!/usr/bin/env python3
"""Generate MkDocs example pages from the repository EXAMPLES/**/README.md.

Usage: python3 tools/docs/build_examples.py --repo . --output build/docs-generated
"""
import argparse
import json
import posixpath
import re
import shutil
import textwrap
from pathlib import Path
from urllib.parse import quote

import yaml

REQUIRED = ("title", "summary", "category", "tags")
SECTIONS = ("Overview", "Main parameters", "How to run", "Expected output")

# Current user-facing explanations supersede the development specifications.
# Example README files retain their original relative links in the repository.
SITE_DOC_LINKS = {
    "docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md": "model/two-layer.md",
    "docs/TWO_LAYER_HYPERBOLICITY.md": "model/hyperbolicity.md",
    "docs/TL13B_C_PHYSICAL_CLOSURES.md": "model/interlayer-processes.md",
    "docs/NETCDF_OUTPUT.md": "getting-started/netcdf-fields.md",
    "docs/RESTART_INPUT.md": "getting-started/initial-conditions.md",
    "docs/website/getting-started/visualization.md": "getting-started/visualization.md",
    "docs/website/getting-started/netcdf-fields.md": "getting-started/netcdf-fields.md",
}



def parse_readme(path):
    """Preserve authored metadata, or derive it for plain Markdown example guides."""
    raw = path.read_text(encoding="utf-8-sig")
    match = re.match(r"\A---\s*\n(.*?)\n---\s*\n", raw, re.S)
    if match:
        meta = yaml.safe_load(match.group(1))
        if not isinstance(meta, dict):
            raise ValueError(f"{path}: invalid metadata")
        for key in REQUIRED:
            if not meta.get(key):
                raise ValueError(f"{path}: missing {key}")
        if not isinstance(meta["tags"], list) or not all(isinstance(t, str) for t in meta["tags"]):
            raise ValueError(f"{path}: tags must be a string list")
        body = raw[match.end():]
    else:
        body = raw
        heading = re.search(r"^#\s+(.+?)\s*$", body, re.M)
        if not heading:
            raise ValueError(f"{path}: missing top-level heading")
        title = heading.group(1).strip()
        paragraphs = re.split(r"\n\s*\n", body[heading.end():])
        description = next(
            (p for p in paragraphs if p.strip() and not p.lstrip().startswith(("#", "|"))),
            title,
        )
        summary = textwrap.shorten(" ".join(description.split()), width=240, placeholder="...")
        folder = path.parent.as_posix().upper()
        if "TWO_LAYER" in folder:
            category = "Two-layer flows"
        elif "EXAMPLE_BUMP" in folder:
            category = "Analytical benchmarks"
        elif any(label in folder for label in ("ETNA", "TAAL", "LAHAR", "VESUVIO")):
            category = "Volcanic flows"
        else:
            category = "Other examples"
        tags = [part.lower() for part in path.parent.name.split("_") if part and part != "EXAMPLE"]
        meta = {"title": title, "summary": summary, "category": category, "tags": tags or ["simulation"]}
    missing = [section for section in SECTIONS
               if not re.search(r"^#{2,4}\s+" + re.escape(section) + r"\s*$", body, re.I | re.M)]
    return meta, body, missing


def is_example_readme(readme, examples_root):
    """Exclude supporting DEM and template guides from the public example catalogue."""
    parts = readme.parent.relative_to(examples_root).parts
    if not parts or not parts[0].startswith("EXAMPLE_"):
        return False
    if len(parts) == 1:
        return True
    if any(part.casefold() in {"dem", "initial_conditions", "templatedir", "output"}
           for part in parts[1:]):
        return False
    return any(
        item.is_file() and item.name.casefold() in {"imex_sflow2d.inp", "imex_sflow2d.template"}
        for item in readme.parent.iterdir()
    )


# Developer-only notes remain in the repository README, not the user manual.
# Keep scientific limitations, run instructions and case bibliographies visible.
INTERNAL_EXAMPLE_SECTIONS = {
    "documentation provenance",
    "documentation verification note",
    "provenance and documentation notes",
    "physical-model namelists in the legacy example",
    "positivity and lie-source diagnostics",
    "l2 face-base and reconstruction audit",
    "wet/dry-front mass diagnostic",
}


def public_example_body(body):
    """Exclude internal migration/diagnostic appendices from website examples."""
    output = []
    hidden = False
    for line in body.splitlines(keepends=True):
        section = re.match(r"^##\s+(.+?)\s*$", line)
        if section:
            hidden = section.group(1).strip().casefold() in INTERNAL_EXAMPLE_SECTIONS
        if not hidden:
            output.append(line)
    return "".join(output).rstrip() + "\n"


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
        if marker != "!" and clean in SITE_DOC_LINKS:
            target = posixpath.relpath(
                SITE_DOC_LINKS[clean], start=f"examples/{example_dir.as_posix()}"
            )
            return f"[{label}]({target})"
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
        if not is_example_readme(readme, examples_root):
            continue
        relative_dir = readme.parent.relative_to(examples_root)
        meta, body, missing = parse_readme(readme)
        body = public_example_body(body)
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
    lines = [
        "# Examples",
        "",
        "Explore simulations organized by application. Each guide describes its "
        "physical setup, required inputs, preparation commands and expected outputs.",
        "",
        "Example parameter values are illustrative unless calibration and "
        "comparison with observations are explicitly documented.",
        "",
    ]
    for category in sorted(categories, key=str.casefold):
        lines.extend([f"## {category}", ""])
        for slug, meta, _ in sorted(categories[category], key=lambda entry: entry[1]["title"].casefold()):
            link = f"{slug}/index.md"
            lines.extend([f"### [{meta['title']}]({link})", "", str(meta["summary"]), "", "Tags: " + ", ".join(f"`{t}`" for t in meta["tags"]), ""])
    (pages / "index.md").write_text("\n".join(lines), encoding="utf-8")
    nav = [
        "  - Home: index.md",
        "  - Getting started:",
        "      - Overview: getting-started/index.md",
        "      - Build from source: getting-started/installation.md",
        "      - Docker: getting-started/docker.md",
        "      - Initial conditions and restart: getting-started/initial-conditions.md",
        "      - Output files: getting-started/output-files.md",
        "      - NetCDF variables and units: getting-started/netcdf-fields.md",
        "      - Visualization (ParaView and QGIS): getting-started/visualization.md",
        "      - Python post-processing: getting-started/postprocessing.md",
        "  - Model:",
        "      - Overview: model/index.md",
        "      - Numerical methods: model/numerics.md",
        "      - Two-layer formulation: model/two-layer.md",
        "      - Coupled pressure: model/coupling.md",
        "      - Interlayer processes: model/interlayer-processes.md",
        "      - Characteristic limits: model/hyperbolicity.md",
        "  - Examples:",
        "      - Catalogue: examples/index.md",
    ]
    for category in sorted(categories, key=str.casefold):
        nav.append(f"      - {json.dumps(category, ensure_ascii=False)}:")
        for slug, meta, _ in sorted(categories[category], key=lambda entry: entry[1]["title"].casefold()):
            label = json.dumps(meta['title'], ensure_ascii=False)
            nav.append(f"          - {label}: examples/{slug}/index.md")
    nav += ["  - References and citation: references/index.md"]
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
