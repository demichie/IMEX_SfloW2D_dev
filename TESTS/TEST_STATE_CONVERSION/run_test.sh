#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-state-conversion.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace -I"$repo_dir/src/state" \
    -c \
    "$repo_dir/src/core/parameters_2d.f90" \
    "$repo_dir/src/core/diagnostics_2d.f90" \
    "$repo_dir/src/core/constitutive_parameters_2d.f90" \
    "$repo_dir/src/numerics/equation_metadata_2d.f90" \
    "$repo_dir/src/numerics/source_catalog_2d.f90" \
    "$repo_dir/src/numerics/source_registry_2d.f90" \
    "$repo_dir/src/core/complexify.f90" \
    "$repo_dir/src/geometry/geometry_2d.f90" \
    "$repo_dir/src/state/state_conversion_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_analytic_utils_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_parameters_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_voellmy_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_plastic_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_froude_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_lucas_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_lava_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_lahar_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_empirical_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_quadratic_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_mui_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_mui_voellmy_2d.f90" \
    "$repo_dir/src/physics/rheology_laws/rheology_2d.f90" \
    "$repo_dir/src/physics/thermal_laws_2d.f90" \
    "$repo_dir/src/physics/pore_pressure_laws_2d.f90" \
    "$repo_dir/src/physics/mass_exchange_laws_2d.f90" \
    "$repo_dir/src/physics/equation_terms_2d.f90" \
    "$test_dir/test_state_conversion.f90"

gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace \
    ./*.o -llapack \
    -o test_state_conversion

./test_state_conversion
