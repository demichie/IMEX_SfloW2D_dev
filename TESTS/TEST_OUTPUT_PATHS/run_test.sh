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

prepare_solver_input() {
    fixture_dir=$1
    if [ -f "$fixture_dir/IMEX_sfloW2D.inp" ]; then
        mv "$fixture_dir/IMEX_sfloW2D.inp" "$fixture_dir/IMEX_input.tmp"
        mv "$fixture_dir/IMEX_input.tmp" "$fixture_dir/IMEX_SfloW2D.inp"
    fi
}

multi_run_dir="$work_dir/multi-run"
mkdir "$multi_run_dir"
cp -R "$test_dir/../TEST_2D/." "$multi_run_dir/"
chmod -R u+w "$multi_run_dir"
prepare_solver_input "$multi_run_dir"
source_input="$multi_run_dir/original.inp"
cp "$multi_run_dir/IMEX_SfloW2D.inp" "$source_input"
for run_name in path-run-a path-run-b; do
    sed "s/RUN_NAME=\"TEST_2D\"/RUN_NAME=\"$run_name\"/" \
        "$source_input" > "$multi_run_dir/input.tmp"
    mv "$multi_run_dir/input.tmp" "$multi_run_dir/IMEX_SfloW2D.inp"
    if ! (
        cd "$multi_run_dir"
        "$build_src_dir/IMEX_SfloW2D" > "$run_name.log" 2>&1
    ); then
        cat "$multi_run_dir/$run_name.log" >&2
        exit 1
    fi
    [ -s "$multi_run_dir/OUTPUT/$run_name/metadata/output_manifest.csv" ] || {
        echo "FAIL: output manifest is missing for $run_name" >&2
        cat "$multi_run_dir/$run_name.log" >&2
        find "$multi_run_dir" -maxdepth 5 -print >&2
        exit 1
    }
done
[ -d "$multi_run_dir/OUTPUT/path-run-a" ] && \
    [ -d "$multi_run_dir/OUTPUT/path-run-b" ] || {
    echo "FAIL: distinct RUN_NAME values did not produce isolated output trees" >&2
    exit 1
}
[ "$(find "$multi_run_dir/OUTPUT" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')" = 2 ] || {
    echo "FAIL: unexpected run directory in shared OUTPUT_BASE_DIRECTORY" >&2
    exit 1
}
echo "PASS: distinct runs remain isolated under a shared output base"
