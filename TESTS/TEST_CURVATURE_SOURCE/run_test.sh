#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-curvature-source.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -I"$repo_dir/src/state" -c \
    "$repo_dir/src/core/parameters_2d.f90" \
    "$repo_dir/src/core/diagnostics_2d.f90" \
    "$repo_dir/src/numerics/lapack_interfaces_2d.f90" \
    "$repo_dir/src/core/constitutive_parameters_2d.f90" \
    "$repo_dir/src/numerics/equation_metadata_2d.f90" \
    "$repo_dir/src/core/complexify.f90" \
    "$repo_dir/src/geometry/geometry_2d.f90" \
    "$repo_dir/src/state/state_conversion_2d.f90" \
    "$repo_dir/src/physics/equation_terms_2d.f90" \
    "$test_dir/test_curvature_source.f90"

gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow \
    parameters_2d.o diagnostics_2d.o lapack_interfaces_2d.o \
    constitutive_parameters_2d.o equation_metadata_2d.o complexify.o \
    geometry_2d.o state_conversion_2d.o equation_terms_2d.o \
    test_curvature_source.o -llapack -o test_curvature_source

./test_curvature_source
