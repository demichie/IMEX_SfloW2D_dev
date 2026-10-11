#!/usr/bin/env python3
"""Prepare the original Etna DEM without avalanche preprocessing.

The archive is shared with EXAMPLE_ETNA. Extraction preserves its exact bytes;
no interpolation, smoothing, ellipsoid excavation or initial-volume generation
is performed here. The solver retains its own normal topography import.
Existing input/DEM files are never overwritten. Only Python's standard library
is required. --directory initializes an isolated execution directory for tests.
"""

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import zipfile


EXAMPLE = Path(__file__).resolve().parent
ARCHIVE = EXAMPLE.parent / "EXAMPLE_ETNA/DEM/Etna2014_crop.zip"
ARCHIVE_SHA256 = "b2fc0de68cede0f432b7cd31c940499fc96515f9716ee49d6df7653711392277"
DEM_SHA256 = "072cfd4961f08bf858afc25c76fc5d37c19e52da477ec0cc95392ab881587bae"
MEMBER = "Etna2014_crop.asc"


def digest(path):
    """Return a streaming SHA-256 without loading the 102 MB raster into memory."""
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def install_dem(archive, destination, archive_hash, member_hash):
    """Install the named raw member atomically, accepting only identical existing data.

    archive/destination are explicit paths; hashes identify the required archive
    and unprocessed raster. Failures raise an exception without changing an
    existing destination. Named-member streaming avoids ZIP path traversal.
    """
    destination = Path(destination)
    if digest(archive) != archive_hash:
        raise ValueError("DEM archive checksum does not match the original Etna fixture")
    if destination.is_symlink():
        raise ValueError("Refusing a symbolic-link DEM destination")
    if destination.exists():
        if not destination.is_file() or digest(destination) != member_hash:
            raise ValueError("Existing DEM differs; it will not be overwritten")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with zipfile.ZipFile(archive) as source, source.open(MEMBER) as stream:
            with tempfile.NamedTemporaryFile(dir=destination.parent,
                                             prefix=".original-dem-", delete=False) as output:
                temporary = Path(output.name)
                shutil.copyfileobj(stream, output, length=1024 * 1024)
        if digest(temporary) != member_hash:
            raise ValueError("Extracted DEM checksum differs from the original raster")
        # Unlike replace(), link() cannot overwrite a concurrently created file.
        os.link(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def prepare(directory):
    """Initialize input and raw DEM in directory, preserving all supplied parameters."""
    directory = Path(directory).resolve()
    original = EXAMPLE / "IMEX_SfloW2D.inp"
    input_path = directory / original.name
    directory.mkdir(parents=True, exist_ok=True)
    if input_path.is_symlink():
        raise ValueError("Refusing a symbolic-link input destination")
    if input_path.exists():
        if not input_path.is_file() or digest(input_path) != digest(original):
            raise ValueError("Existing input differs; it will not be overwritten")
    else:
        with original.open("rb") as source, input_path.open("xb") as output:
            shutil.copyfileobj(source, output)
    install_dem(ARCHIVE, directory / "DEM" / MEMBER, ARCHIVE_SHA256, DEM_SHA256)
    return directory / "DEM" / MEMBER


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=EXAMPLE,
                        help="execution directory; default is this example directory")
    print("Verified original DEM:", prepare(parser.parse_args().directory))
