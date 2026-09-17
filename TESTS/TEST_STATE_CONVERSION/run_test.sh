#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-state-conversion.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -fcheck=all -fbacktrace \
    -c \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/equation_metadata_2d.f90" \
    "$repo_dir/src/complexify.f90" \
    "$repo_dir/src/geometry_2d.f90" \
    "$repo_dir/src/constitutive_2d.f90" \
    "$test_dir/test_state_conversion.f90"

gfortran -O0 -g -fcheck=all -fbacktrace \
    parameters_2d.o equation_metadata_2d.o complexify.o geometry_2d.o \
    constitutive_2d.o test_state_conversion.o -llapack \
    -o test_state_conversion

./test_state_conversion
