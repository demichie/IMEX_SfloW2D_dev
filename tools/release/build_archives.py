#!/usr/bin/env python3
"""Build versioned source and Examples archives from one Git tag."""

import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path, PurePosixPath


VERSION_PATTERN = re.compile(r"v?\d+(?:\.\d+)*(?:[-+][0-9A-Za-z.-]+)?\Z")
AC_INIT_PATTERN = re.compile(
    r"(?m)^(AC_INIT\(\[IMEX_SfloW2D\], )\[[^\]]+\](?=,)"
)
AC_INIT_VERSION_PATTERN = re.compile(
    r"(?m)^AC_INIT\(\[IMEX_SfloW2D\], \[([^\]]+)\]"
)
CONFIGURE_VERSION_PATTERN = re.compile(r"(?m)^PACKAGE_VERSION='([^']+)'$")
# Keep this list aligned with Makefile.am's EXTRA_DIST entries for make check.
MAKE_CHECK_FILES = (
    "TESTS/run_tests.sh",
    "TESTS/TEST_OUTPUT_PATHS/run_test.sh",
    "TESTS/TEST_OUTPUT_PATHS/test_snapshot_index.f90",
    "TESTS/TEST_2D/IMEX_sfloW2D.inp",
    "TESTS/TEST_2D/INIT_2D_0000.q_2d",
    "TESTS/TEST_2D/cleanFolder.sh",
    "TESTS/TEST_2D/topography_dem.asc",
    "TESTS/TEST_BOX/IMEX_sfloW2D.inp",
    "TESTS/TEST_BOX/INIT_BOX_0000.q_2d",
    "TESTS/TEST_BOX/cleanFolder.sh",
    "TESTS/TEST_BOX/topography_dem.asc",
    "TESTS/TEST_BW/IMEX_sfloW2D.inp",
    "TESTS/TEST_BW/INIT_BW_0000.q_2d",
    "TESTS/TEST_BW/cleanFolder.sh",
    "TESTS/TEST_BW/topography_dem.asc",
    "TESTS/TEST_CF/IMEX_sfloW2D.inp",
    "TESTS/TEST_CF/cleanFolder.sh",
    "TESTS/TEST_TAAL/IMEX_sfloW2D.inp",
    "TESTS/TEST_TAAL/cleanFolder.sh",
    "TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/IMEX_SfloW2D.inp",
    "TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/cleanFolder.sh",
    "TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/example_2D_0000.q_2d",
    "TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/topography_dem.asc",
    "TESTS/TEST_DAMBREAK/TEST_OFF_FINH/IMEX_SfloW2D.inp",
    "TESTS/TEST_DAMBREAK/TEST_OFF_FINH/cleanFolder.sh",
    "TESTS/TEST_DAMBREAK/TEST_OFF_FINH/example_2D_0000.q_2d",
    "TESTS/TEST_DAMBREAK/TEST_OFF_FINH/topography_dem.asc",
    "TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/IMEX_SfloW2D.inp",
    "TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/cleanFolder.sh",
    "TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/example_2D_0000.q_2d",
    "TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/topography_dem.asc",
    "TESTS/TEST_DAMBREAK/TEST_STEP_FINH/IMEX_SfloW2D.inp",
    "TESTS/TEST_DAMBREAK/TEST_STEP_FINH/cleanFolder.sh",
    "TESTS/TEST_DAMBREAK/TEST_STEP_FINH/example_2D_0000.q_2d",
    "TESTS/TEST_DAMBREAK/TEST_STEP_FINH/topography_dem.asc",
    "EXAMPLES/EXAMPLE_CF/topo.zip",
    "EXAMPLES/EXAMPLE_TAAL/topography_dem.asc",
)


def run(command, cwd=None, capture=False):
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
    )


def git_output(repo, *args):
    return run(["git", *args], cwd=repo, capture=True).stdout.strip()


def validate_members(members, archive_name):
    for member in members:
        name = PurePosixPath(member.name)
        if name.is_absolute() or ".." in name.parts:
            raise ValueError(f"Unsafe archive path in {archive_name}: {member.name}")
        if member.issym() or member.islnk():
            link_target = PurePosixPath(member.linkname)
            if link_target.is_absolute():
                raise ValueError(f"Absolute link in {archive_name}: {member.name}")
            combined = name.parent / link_target
            depth = 0
            for part in combined.parts:
                if part == "..":
                    depth -= 1
                elif part not in ("", "."):
                    depth += 1
                if depth < 0:
                    raise ValueError(f"Escaping link in {archive_name}: {member.name}")
    return members


def archive_members(path):
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
    return validate_members(members, path.name)


def safe_extract(path, destination):
    members = archive_members(path)
    with tarfile.open(path, "r:gz") as archive:
        archive.extractall(destination, members=members)


def tracked_paths(repo, revision, pathspec):
    output = git_output(repo, "ls-tree", "-r", "--name-only", revision, "--", pathspec)
    return set(output.splitlines()) if output else set()


def assert_source_inventory(repo, revision, source_archive, version):
    members = archive_members(source_archive)
    paths = {member.name for member in members}
    prefix = f"IMEX_SfloW2D-{version}/"
    required = {
        "configure",
        "configure.ac",
        "Makefile.in",
        "src/Makefile.in",
        "COPYING",
        "README.md",
        *MAKE_CHECK_FILES,
    }
    missing = sorted(prefix + name for name in required if prefix + name not in paths)
    if missing:
        raise ValueError("Source archive is missing required files: " + ", ".join(missing))

    source_files = tracked_paths(repo, revision, "src")
    required_source = {
        path for path in source_files
        if path.endswith((".f90", ".F90", ".inc", ".c", ".h"))
    }
    missing_source = sorted(prefix + name for name in required_source if prefix + name not in paths)
    if missing_source:
        raise ValueError("Source archive is missing source files: " + ", ".join(missing_source))
    with tarfile.open(source_archive, "r:gz") as archive:
        configure_ac = archive.extractfile(prefix + "configure.ac")
        configure = archive.extractfile(prefix + "configure")
        if configure_ac is None or configure is None:
            raise ValueError("Source archive is missing Autotools version metadata")
        source_version = AC_INIT_VERSION_PATTERN.search(
            configure_ac.read().decode("utf-8")
        )
        generated_version = CONFIGURE_VERSION_PATTERN.search(
            configure.read().decode("utf-8")
        )
    if not source_version or source_version.group(1) != version:
        raise ValueError("configure.ac package version does not match the tag")
    if not generated_version or generated_version.group(1) != version:
        raise ValueError("generated configure package version does not match the tag")
    return sorted(member.name for member in members)


def assert_examples_inventory(repo, revision, examples_archive, version, extract_dir):
    members = archive_members(examples_archive)
    prefix = f"IMEX_SfloW2D-{version}-examples/"
    tracked = tracked_paths(repo, revision, "EXAMPLES")
    expected = {prefix + name for name in tracked}
    actual_files = {member.name for member in members if not member.isdir()}
    if expected != actual_files:
        missing = sorted(expected - actual_files)
        unexpected = sorted(actual_files - expected)
        raise ValueError(
            f"Examples archive inventory differs from tag; missing={missing}, "
            f"unexpected={unexpected}"
        )

    readmes = sorted(path for path in tracked if path.endswith("/README.md"))
    if not readmes:
        raise ValueError("No standardized EXAMPLES/**/README.md files found at tag")
    safe_extract(examples_archive, extract_dir)
    missing_readmes = [
        path for path in readmes if not (extract_dir / prefix / path).is_file()
    ]
    if missing_readmes:
        raise ValueError("Examples extraction is missing README files: " + ", ".join(missing_readmes))
    return sorted(member.name for member in members), readmes


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True, help="release tag to package")
    parser.add_argument("--output", required=True, type=Path, help="new or empty output directory")
    parser.add_argument(
        "--configure-option",
        action="append",
        default=[],
        help="option passed to configure; may be repeated",
    )
    args = parser.parse_args()
    if not VERSION_PATTERN.fullmatch(args.tag):
        parser.error("--tag must be a version-like tag, such as 1.2.3 or v1.2.3")

    repo = Path(git_output(Path.cwd(), "rev-parse", "--show-toplevel"))
    revision = git_output(repo, "rev-parse", "--verify", f"refs/tags/{args.tag}^{{commit}}")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_name = f"IMEX_SfloW2D-{args.tag}.tar.gz"
    examples_name = f"IMEX_SfloW2D-{args.tag}-examples.tar.gz"
    manifest_name = f"IMEX_SfloW2D-{args.tag}-manifest.json"
    destinations = [output / name for name in (source_name, examples_name, manifest_name)]
    existing = [path.name for path in destinations if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite: " + ", ".join(existing))

    with tempfile.TemporaryDirectory(prefix="imex-release-") as temporary:
        temporary = Path(temporary)
        staged = [temporary / name for name in (source_name, examples_name, manifest_name)]
        source_tar = temporary / "source.tar"
        run(
            [
                "git", "archive", "--format=tar",
                f"--prefix=IMEX_SfloW2D-{args.tag}/",
                f"--output={source_tar}", revision,
            ],
            cwd=repo,
        )
        source_tree = temporary / f"IMEX_SfloW2D-{args.tag}"
        with tarfile.open(source_tar, "r:") as archive:
            members = validate_members(archive.getmembers(), source_tar.name)
            archive.extractall(temporary, members=members)

        configure_ac = source_tree / "configure.ac"
        configure_text = configure_ac.read_text(encoding="utf-8")
        configure_text, replaced = AC_INIT_PATTERN.subn(
            rf"\g<1>[{args.tag}]",
            configure_text,
        )
        if replaced != 1:
            raise ValueError("Could not set the Autotools package version from the release tag")
        configure_ac.write_text(configure_text, encoding="utf-8")
        run(["autoreconf", "-fi"], cwd=source_tree)

        build = temporary / "build"
        build.mkdir()
        run([str(source_tree / "configure"), *args.configure_option], cwd=build)
        run(
            [
                "make",
                "dist",
                f"VERSION={args.tag}",
                "PACKAGE=IMEX_SfloW2D",
                "PACKAGE_TARNAME=IMEX_SfloW2D",
            ],
            cwd=build,
        )
        built_source = build / source_name
        if not built_source.is_file():
            raise FileNotFoundError(f"Autotools did not create {built_source}")
        shutil.copy2(built_source, staged[0])

        examples_tar = temporary / "examples.tar"
        run(
            [
                "git", "archive", "--format=tar",
                f"--prefix=IMEX_SfloW2D-{args.tag}-examples/",
                f"--output={examples_tar}", revision, "EXAMPLES",
            ],
            cwd=repo,
        )
        with examples_tar.open("rb") as source, staged[1].open("wb") as target:
            with gzip.GzipFile(filename="", mode="wb", fileobj=target, mtime=0) as compressed:
                shutil.copyfileobj(source, compressed)

        source_inventory = assert_source_inventory(
            repo, revision, staged[0], args.tag
        )
        extracted_examples = temporary / "extracted-examples"
        extracted_examples.mkdir()
        examples_inventory, readmes = assert_examples_inventory(
            repo, revision, staged[1], args.tag, extracted_examples
        )

        manifest = {
            "schema_version": 1,
            "tag": args.tag,
            "package_version": args.tag,
            "git_revision": revision,
            "archives": {
                "source": {
                    "filename": source_name,
                    "sha256": sha256(staged[0]),
                    "inventory": source_inventory,
                },
                "examples": {
                    "filename": examples_name,
                    "sha256": sha256(staged[1]),
                    "inventory": examples_inventory,
                    "standardized_readmes": readmes,
                },
            },
        }
        staged[2].write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        for source_path, destination in zip(staged, destinations):
            shutil.copy2(source_path, destination)
    print(f"Built {source_name}")
    print(f"Built {examples_name}")
    print(f"Built {manifest_name}")
    print(f"Revision: {revision}")
    print(f"Standardized example READMEs: {len(readmes)}")


if __name__ == "__main__":
    main()
