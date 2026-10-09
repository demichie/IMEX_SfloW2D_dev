#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-rheology-froude.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

# Trap invalid arithmetic in both profiles: optimized success alone must not
# conceal a zero/negative Froude denominator. No scientific tolerance is relaxed.
for profile in strict optimized; do
    case "$profile" in
        strict) profile_flags='-O0 -g -fcheck=all' ;;
        optimized) profile_flags='-Ofast -funroll-all-loops' ;;
    esac
    mkdir "$profile"
    cd "$profile"
    gfortran $profile_flags -fopenmp -fbacktrace \
        -ffpe-trap=invalid,zero,overflow -c \
        "$repo_dir/src/parameters_2d.f90" \
        "$repo_dir/src/diagnostics_2d.f90" \
        "$repo_dir/src/constitutive_parameters_2d.f90" \
        "$repo_dir/src/equation_metadata_2d.f90" \
        "$repo_dir/src/complexify.f90" \
        "$repo_dir/src/geometry_2d.f90" \
        "$repo_dir/src/state_conversion_2d.f90" \
        "$repo_dir/src/equation_terms_2d.f90"

    # Keep assertion/IEEE checks strict even when the numerical modules use
    # -Ofast, which otherwise permits optimizing away NaN/Inf checks.
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace \
        -ffpe-trap=invalid,zero,overflow -c "$test_dir/test_rheology_froude.f90"

    gfortran $profile_flags -fopenmp -fbacktrace \
        -ffpe-trap=invalid,zero,overflow \
        parameters_2d.o diagnostics_2d.o constitutive_parameters_2d.o \
        equation_metadata_2d.o complexify.o geometry_2d.o \
        state_conversion_2d.o equation_terms_2d.o test_rheology_froude.o \
        -llapack -o test_rheology_froude

    ./test_rheology_froude
    cd ..
done
