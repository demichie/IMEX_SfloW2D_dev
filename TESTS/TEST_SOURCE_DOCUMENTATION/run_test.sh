#!/bin/sh
# Validate documentation coverage without a compiler or a Doxygen installation.
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/check_documentation.py" "$@"
