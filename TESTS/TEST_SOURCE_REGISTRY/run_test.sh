#!/bin/sh
# D17-B: metadata-only checks; all build products stay outside src.
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
work_dir=${D17_WORK_DIR:-$(mktemp -d /tmp/imex-d17-registry.XXXXXX)}
mkdir -p "$work_dir"
work_dir=$(CDPATH= cd -- "$work_dir" && pwd)
compiler=${FC:-gfortran}
for profile in strict optimized; do
    mkdir -p "$work_dir/$profile"
    cd "$work_dir/$profile"
    case "$profile" in
      strict) flags='-O0 -g -Wall -Wextra -fcheck=all -finit-real=snan' ;;
      optimized) flags='-Ofast -funroll-all-loops' ;;
    esac
    "$compiler" $flags -fbacktrace -ffpe-trap=invalid,zero,overflow -c \
      "$repo_dir/src/core/parameters_2d.f90" \
      "$repo_dir/src/core/constitutive_parameters_2d.f90" \
      "$repo_dir/src/numerics/equation_metadata_2d.f90" \
      "$repo_dir/src/numerics/source_catalog_2d.f90" \
      "$repo_dir/src/numerics/source_input_2d.f90" \
      "$repo_dir/src/numerics/source_registry_2d.f90" \
      "$repo_dir/src/physics/source_configuration_2d.f90" \
      "$test_dir/test_source_registry.f90"
    "$compiler" -fbacktrace -ffpe-trap=invalid,zero,overflow ./*.o -o test_source_registry
    ./test_source_registry > results.log
    cat results.log
done
printf 'D17 evidence directory: %s\n' "$work_dir"
