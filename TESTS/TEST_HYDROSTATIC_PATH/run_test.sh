#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-hydrostatic-path.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -fcheck=all -fbacktrace -c \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/geometry_2d.f90" \
    "$repo_dir/src/well_balanced_2d.f90" \
    "$test_dir/test_hydrostatic_path.f90"

gfortran -O0 -g -fcheck=all -fbacktrace \
    parameters_2d.o geometry_2d.o well_balanced_2d.o \
    test_hydrostatic_path.o -llapack -o test_hydrostatic_path

./test_hydrostatic_path
