#!/bin/sh
# Validate documentation coverage without a compiler or a Doxygen installation.
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
PYTHONDONTWRITEBYTECODE=1 python3 "$repo_dir/tools/docs/check_fortran_style.py" "$repo_dir/src"
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/check_documentation.py" "$@"
