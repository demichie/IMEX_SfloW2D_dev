#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
python3 -B "$test_dir/test_projection_checker.py"
build_dir=$(mktemp -d /tmp/imex-mass-exchange.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $build_dir"
else
    trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
fi

netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)

# Compile the production module order, excluding only the application main.
# Real checkpoint routines are exercised, not a test-only serialization stub.
sources=$(awk '
    /^IMEX_SfloW2D_SOURCES =/ { in_sources=1; next }
    in_sources && /^[[:space:]]*$/ { exit }
    in_sources { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
' "$repo_dir/src/Makefile.am")
production_objects=''
for source in $sources; do
    production_objects="$production_objects ${source%.f90}.o"
done

for profile in strict optimized; do
    mkdir "$build_dir/$profile"
    cd "$build_dir/$profile"
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all' ;;
        optimized) flags='-Ofast -funroll-all-loops' ;;
    esac
    for source in $sources; do
        gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            $netcdf_fflags -c "$repo_dir/src/$source"
    done
    # Assertion/reference arithmetic must not inherit finite-math assumptions.
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace \
        -ffpe-trap=invalid,zero,overflow $netcdf_fflags \
        -c "$test_dir/test_mass_exchange.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        ./*.o -llapack $netcdf_flibs $netcdf_clibs -o test_mass_exchange
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace \
        -ffpe-trap=invalid,zero,overflow $netcdf_fflags \
        -c "$test_dir/test_projection_stress.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $production_objects test_projection_stress.o \
        -llapack $netcdf_flibs $netcdf_clibs -o test_projection_stress

    for threads in 1 4; do
        mkdir "threads-$threads"
        (
            cd "threads-$threads"
            N8_CASE=all OMP_NUM_THREADS=$threads OMP_DYNAMIC=FALSE \
                ../test_mass_exchange "$threads" > test.log
            python3 "$test_dir/check_diagnostics.py" test.log "$profile" "$threads"
            OMP_NUM_THREADS=$threads OMP_DYNAMIC=FALSE \
                ../test_projection_stress "$threads" > projection.log
            python3 "$test_dir/check_projection_stress.py" projection.log "$profile" "$threads" \
                > projection_evidence.log
        )
    done
    cmp threads-1/snapshot.bin threads-4/snapshot.bin
    cmp threads-1/checkpoint.bin threads-4/checkpoint.bin
    cmp threads-1/gas_snapshot.bin threads-4/gas_snapshot.bin
    cmp threads-1/flat_evaluations.bin threads-4/flat_evaluations.bin
    cmp threads-1/projection_snapshot.bin threads-4/projection_snapshot.bin
    echo "PASS: $profile mass exchange and checkpoint continuation, actual 1/4 threads"
done
