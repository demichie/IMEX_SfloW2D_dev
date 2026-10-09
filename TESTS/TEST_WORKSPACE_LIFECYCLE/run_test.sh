#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-workspace-lifecycle.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
cd "$build_dir"

gfortran -O0 -g -Wall -Wextra -fopenmp -fcheck=all -fbacktrace -c \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/diagnostics_2d.f90" \
    "$repo_dir/src/constitutive_parameters_2d.f90" \
    "$repo_dir/src/model_layout_2d.f90" \
    "$repo_dir/src/equation_metadata_2d.f90" \
    "$repo_dir/src/complexify.f90" \
    "$repo_dir/src/geometry_2d.f90" \
    "$repo_dir/src/state_conversion_2d.f90" \
    "$repo_dir/src/equation_terms_2d.f90" \
    "$repo_dir/src/domain_2d.f90" \
    "$repo_dir/src/state_2d.f90" \
    "$repo_dir/src/hp_reconstruction_2d.f90" \
    "$repo_dir/src/reconstruction_2d.f90" \
    "$repo_dir/src/nonconservative_2d.f90" \
    "$repo_dir/src/pccu_2d.f90" \
    "$repo_dir/src/hyperbolic_2d.f90" \
    "$repo_dir/src/mass_exchange_2d.f90" \
    "$repo_dir/src/init_2d.f90" \
    "$test_dir/test_workspace_lifecycle.f90"
gfortran -fopenmp -fcheck=all -fbacktrace ./*.o -llapack -o test_workspace_lifecycle
OMP_NUM_THREADS=4 ./test_workspace_lifecycle
