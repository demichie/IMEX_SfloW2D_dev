#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-workspace-lifecycle.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
cd "$build_dir"

gfortran -O0 -g -Wall -Wextra -fopenmp -fcheck=all -fbacktrace -I"$repo_dir/src/state" -c \
    "$repo_dir/src/core/parameters_2d.f90" \
    "$repo_dir/src/core/diagnostics_2d.f90" \
    "$repo_dir/src/core/constitutive_parameters_2d.f90" \
    "$repo_dir/src/core/model_layout_2d.f90" \
    "$repo_dir/src/numerics/equation_metadata_2d.f90" \
    "$repo_dir/src/core/complexify.f90" \
    "$repo_dir/src/geometry/geometry_2d.f90" \
    "$repo_dir/src/state/state_conversion_2d.f90" \
    "$repo_dir/src/physics/equation_terms_2d.f90" \
    "$repo_dir/src/geometry/domain_2d.f90" \
    "$repo_dir/src/state/state_2d.f90" \
    "$repo_dir/src/spatial/hp_reconstruction_2d.f90" \
    "$repo_dir/src/spatial/reconstruction_2d.f90" \
    "$repo_dir/src/spatial/nonconservative_2d.f90" \
    "$repo_dir/src/spatial/pccu_2d.f90" \
    "$repo_dir/src/spatial/hyperbolic_2d.f90" \
    "$repo_dir/src/physics/mass_exchange_2d.f90" \
    "$repo_dir/src/core/init_2d.f90" \
    "$test_dir/test_workspace_lifecycle.f90"
gfortran -fopenmp -fcheck=all -fbacktrace ./*.o -llapack -o test_workspace_lifecycle
OMP_NUM_THREADS=4 ./test_workspace_lifecycle
