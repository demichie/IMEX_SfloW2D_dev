#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-state-initialization.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -Wall -Wextra -fopenmp -fcheck=all -fbacktrace \
    -c \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/diagnostics_2d.f90" \
    "$repo_dir/src/constitutive_parameters_2d.f90" \
    "$repo_dir/src/equation_metadata_2d.f90" \
    "$repo_dir/src/complexify.f90" \
    "$repo_dir/src/geometry_2d.f90" \
    "$repo_dir/src/state_conversion_2d.f90" \
    "$repo_dir/src/equation_terms_2d.f90" \
    "$repo_dir/src/constitutive_2d.f90" \
    "$repo_dir/src/state_2d.f90" \
    "$repo_dir/src/domain_2d.f90" \
    "$test_dir/test_state_initialization.f90"

gfortran -O0 -g -Wall -Wextra -fopenmp -fcheck=all -fbacktrace \
    parameters_2d.o diagnostics_2d.o constitutive_parameters_2d.o \
    equation_metadata_2d.o \
    complexify.o geometry_2d.o state_conversion_2d.o equation_terms_2d.o \
    constitutive_2d.o state_2d.o domain_2d.o test_state_initialization.o \
    -llapack -o test_state_initialization

./test_state_initialization
