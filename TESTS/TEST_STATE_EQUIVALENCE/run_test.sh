#!/bin/sh
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-test-state-equivalence.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
cd "$build_dir"
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)
for source in core/parameters_2d.f90 core/constitutive_parameters_2d.f90 core/diagnostics_2d.f90 state/state_conversion_2d.f90; do
    gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace -I"$repo_dir/src/state" $netcdf_fflags -c "$repo_dir/src/$source"
done
for source in state_types.f90 state_composition.f90 state_mixture.f90 state_conversions.f90 state_physical.f90 state.f90; do
    gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace -I"$repo_dir/src/state" -I"$repo_dir/src/state" $netcdf_fflags -c "$repo_dir/src/state/$source"
done
gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace -I"$repo_dir/src/state" -I"$repo_dir/src/state" $netcdf_fflags -c "$repo_dir/src/physics/equation_terms_2d.f90"
gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace -I"$repo_dir/src/state" -I"$repo_dir/src/state" $netcdf_fflags -c "$test_dir/test_state_equivalence.f90"
gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace ./*.o -llapack -lblas $netcdf_flibs $netcdf_clibs -lm -o test_state_equivalence
./test_state_equivalence
