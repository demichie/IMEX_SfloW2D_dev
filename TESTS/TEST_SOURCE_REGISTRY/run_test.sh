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
    # Follow production compilation order: source_input_2d depends on
    # physical mass-exchange modules, including lava deposition.
    sources=$(awk '
      /^IMEX_SfloW2D_SOURCES =/ { inside=1; next }
      inside && /^[[:space:]]*$/ { exit }
      inside { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
    ' "$repo_dir/src/Makefile.am")
    netcdf_fflags=$(nf-config --fflags)
    netcdf_flibs=$(nf-config --flibs)
    netcdf_clibs=$(nc-config --libs)
    for source in $sources; do
      case "$source" in
        *.c) cc -c "$repo_dir/src/$source" ;;
        *) "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
             -I"$repo_dir/src/state" $netcdf_fflags -c "$repo_dir/src/$source" ;;
      esac
    done
    "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $netcdf_fflags -c "$test_dir/test_source_registry.f90"
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      ./*.o -llapack $netcdf_flibs $netcdf_clibs -o test_source_registry
    # D17.20 independent Newton-only extension fixture (no historical solver API).
    "$compiler" $flags -fbacktrace -ffpe-trap=invalid,zero,overflow -c \
      "$repo_dir/src/core/diagnostics_2d.f90" \
      "$repo_dir/src/physics/rheology_laws/rheology_newmodel_2d.f90" \
      "$test_dir/test_newmodel_numerics.f90"
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      parameters_2d.o diagnostics_2d.o source_catalog_2d.o \
      rheology_newmodel_2d.o test_newmodel_numerics.o -o test_newmodel_numerics
    ./test_newmodel_numerics > newmodel-numerics.log
    cat newmodel-numerics.log
    ./test_source_registry > results.log
    cat results.log
done
printf 'D17 evidence directory: %s\n' "$work_dir"
