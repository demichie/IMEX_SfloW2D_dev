#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-diagnostics.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -Wall -Wextra -fopenmp -fcheck=all -fbacktrace \
    "$repo_dir/src/diagnostics_2d.f90" \
    "$test_dir/test_diagnostics.f90" \
    -o test_diagnostics

./test_diagnostics noninteractive
printf '\n' | ./test_diagnostics interactive
OMP_NUM_THREADS=2 ./test_diagnostics parallel </dev/null

if ./test_diagnostics fatal >fatal.log 2>&1; then
    echo 'FAIL: fatal_error returned successfully' >&2
    exit 1
fi

grep 'FATAL ERROR: expected test failure' fatal.log >/dev/null
echo 'PASS: fatal diagnostics terminate with failure status'
