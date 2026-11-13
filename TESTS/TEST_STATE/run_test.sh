#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-test-state-new.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
cd "$build_dir"

common="-O0 -g -Wall -Wextra -fcheck=all -fbacktrace -I$repo_dir/src/state"
gfortran $common -c "$repo_dir/src/core/parameters_2d.f90"
gfortran $common -c "$repo_dir/src/core/constitutive_parameters_2d.f90"
gfortran $common -c "$repo_dir/src/core/solid_transport_parameters_2d.f90"
gfortran $common -J. -c "$repo_dir/src/state/state_types.f90"
gfortran $common -J. -c "$repo_dir/src/state/state_composition.f90"
gfortran $common -J. -c "$repo_dir/src/state/state_mixture.f90"
gfortran $common -J. -c "$repo_dir/src/state/state_conversions.f90"
gfortran $common -J. -c "$repo_dir/src/state/state_physical.f90"
gfortran $common -J. -c "$test_dir/test_state.f90"
gfortran $common ./*.o -o test_state
./test_state
