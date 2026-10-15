#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 fortran-compiler openmp-flags /path/to/build/src" >&2
    exit 2
fi

fortran_compiler=$1
openmp_flags=$2
build_src_dir=$3
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/imex-output-path-tests.XXXXXX")
work_dir=$(CDPATH= cd -- "$work_dir" && pwd -P)
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM

# OPENMP_FCFLAGS intentionally supports multiple compiler arguments.
# shellcheck disable=SC2086
"$fortran_compiler" $openmp_flags \
    -I"$build_src_dir" \
    "$test_dir/test_snapshot_index.f90" \
    "$build_src_dir/io/output_paths.o" \
    "$build_src_dir/io/output_dirs.o" \
    "$build_src_dir/core/parameters_2d.o" \
    "$build_src_dir/core/diagnostics_2d.o" \
    -o "$work_dir/test_snapshot_index"

"$work_dir/test_snapshot_index" "$work_dir"
echo "PASS: snapshot index labels include 0, 9999, 10000, and 10001"
