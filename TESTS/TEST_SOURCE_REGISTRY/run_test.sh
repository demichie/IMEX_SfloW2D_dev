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
             -I"$repo_dir/src/state" -I"$repo_dir/src/physics/rheology_laws" -I"$repo_dir/src/physics/thermal_laws" $netcdf_fflags -c "$repo_dir/src/$source" ;;
      esac
    done
    "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $netcdf_fflags -c "$test_dir/test_source_registry.f90"
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      ./*.o -llapack $netcdf_flibs $netcdf_clibs -o test_source_registry
    # The production lava viscous law has an exact linear subcase (zero yield).
    # It keeps the independent four-case Newton/complex-step regression without
    # registering a dummy constitutive law in the solver.
    "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $netcdf_fflags -c "$test_dir/test_lava_newton_numerics.f90"
    production_objects=$(printf '%s\n' $sources | sed -n \
      -e 's|.*/||;s/\.f90$/.o/p' -e 's|.*/||;s/\.c$/.o/p' | tr '\n' ' ')
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $production_objects test_lava_newton_numerics.o -llapack \
      $netcdf_flibs $netcdf_clibs -o test_lava_newton_numerics
    ./test_lava_newton_numerics > lava-newton-numerics.log
    cat lava-newton-numerics.log
    ./test_source_registry > results.log
    cat results.log
    "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $netcdf_fflags -c "$test_dir/test_interlayer_mass_exchange.f90"
    # The paired physical coordinator uses the existing deposition and state
    # modules: link the entire production module closure, excluding test mains.
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $production_objects test_interlayer_mass_exchange.o -llapack \
      $netcdf_flibs $netcdf_clibs -o test_interlayer_mass_exchange
    ./test_interlayer_mass_exchange > interlayer-results.log
    cat interlayer-results.log
    "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $netcdf_fflags -c "$test_dir/test_interlayer_lie_transaction.f90"
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $production_objects test_interlayer_lie_transaction.o -llapack \
      $netcdf_flibs $netcdf_clibs -o test_interlayer_lie_transaction
    if ! ./test_interlayer_lie_transaction > interlayer-lie-results.log 2>&1; then
      printf 'FAIL: interlayer Lie transaction (%s)\n' "$profile" >&2
      cat interlayer-lie-results.log >&2
      exit 1
    fi
    cat interlayer-lie-results.log
    "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $netcdf_fflags -c "$test_dir/test_interlayer_momentum_exchange.f90"
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $production_objects test_interlayer_momentum_exchange.o -llapack \
      $netcdf_flibs $netcdf_clibs -o test_interlayer_momentum_exchange
    ./test_interlayer_momentum_exchange > interlayer-momentum-results.log
    cat interlayer-momentum-results.log
    # Each unused physical method must terminate with the expected fatal
    # diagnostic; merely returning an error status is not sufficient.
    "$compiler" $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $netcdf_fflags -c "$test_dir/test_interlayer_unsupported_treatments.f90"
    "$compiler" -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
      $production_objects test_interlayer_unsupported_treatments.o -llapack \
      $netcdf_flibs $netcdf_clibs -o test_interlayer_unsupported_treatments
    set -- \
      'interlayer_settling: unsupported semi-implicit physical treatment' \
      'interlayer_settling: unsupported implicit physical treatment' \
      'interlayer_settling: unsupported analytic physical treatment' \
      'interlayer_elutriation: unsupported semi-implicit physical treatment' \
      'interlayer_elutriation: unsupported implicit physical treatment' \
      'interlayer_elutriation: unsupported analytic physical treatment' \
      'interlayer_mass_exchange: unsupported semi-implicit physical treatment' \
      'interlayer_mass_exchange: unsupported implicit physical treatment' \
      'interlayer_mass_exchange: unsupported analytic physical treatment' \
      'interlayer_drag_chezy: unsupported explicit physical treatment' \
      'interlayer_drag_chezy: unsupported semi-implicit physical treatment' \
      'interlayer_drag_chezy: unsupported implicit physical treatment' \
      'interlayer_momentum_exchange: unsupported explicit physical treatment' \
      'interlayer_momentum_exchange: unsupported semi-implicit physical treatment' \
      'interlayer_momentum_exchange: unsupported implicit physical treatment'
    case_index=0
    for expected in "$@"; do
      case_index=$((case_index+1))
      if ./test_interlayer_unsupported_treatments "$case_index" > unsupported.log 2>&1; then
        echo "FAIL: unsupported treatment $case_index returned successfully" >&2
        exit 1
      fi
      if ! grep -F "FATAL ERROR:" unsupported.log > /dev/null || \
         ! grep -F "$expected" unsupported.log > /dev/null; then
        cat unsupported.log >&2
        echo "FAIL: unsupported treatment $case_index emitted the wrong diagnostic" >&2
        exit 1
      fi
    done
    printf 'PASS: 15 unsupported interlayer treatments fail with fatal diagnostics (%s)\n' "$profile"
done
printf 'D17 evidence directory: %s\n' "$work_dir"
