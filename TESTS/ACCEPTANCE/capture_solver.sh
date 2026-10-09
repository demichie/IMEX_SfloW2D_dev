#!/bin/sh
# Adapter used only by the acceptance audit inside temporary test fixtures.
set -eu
adapter_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec python3 "$adapter_dir/capture_solver.py" "$@"
