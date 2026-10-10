#!/bin/sh
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/test_checker.py"
build_dir=$(mktemp -d /tmp/imex-n7-geometry-dynamic.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $build_dir"
else
    trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
fi
sources=$(awk '
    /^IMEX_SfloW2D_SOURCES =/ { in_sources=1; next }
    in_sources && /^[[:space:]]*$/ { exit }
    in_sources { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
' "$repo_dir/src/Makefile.am")
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)
failed=0
for profile in strict optimized; do
    mkdir "$build_dir/$profile"
    cd "$build_dir/$profile"
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all' ;;
        optimized) flags='-Ofast -funroll-all-loops' ;;
    esac
    objects=''
    for source in $sources; do
        gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            $netcdf_fflags -c "$repo_dir/src/$source"
        objects="$objects ${source%.f90}.o"
    done
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $netcdf_fflags -c "$repo_dir/TESTS/TEST_IMEX_STAGES/test_imex_stages.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $objects test_imex_stages.o -llapack $netcdf_flibs $netcdf_clibs -o test_imex_stages
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $netcdf_fflags -c "$repo_dir/TESTS/TEST_SPATIAL_OPERATOR/test_spatial_operator.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $objects test_spatial_operator.o -llapack $netcdf_flibs $netcdf_clibs -o test_spatial_operator
    if ! PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/check_cases.py" \
        "$build_dir/$profile/test_imex_stages" "$profile" \
        --spatial-executable "$build_dir/$profile/test_spatial_operator"; then
        failed=1
    fi
done
exit "$failed"
