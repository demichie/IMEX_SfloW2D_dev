#!/bin/sh
# Compile the shared observer against both production profiles; do not change sources.
set -eu
diagnostic=0
equilibrium_only=''
case "${1:-}" in
    '') ;;
    --diagnose) diagnostic=1 ;;
    --equilibrium-only) equilibrium_only='--equilibrium-only' ;;
    *) echo 'Usage: run_test.sh [--diagnose|--equilibrium-only]' >&2; exit 2 ;;
esac
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/test_checker.py"
build_dir=$(mktemp -d /tmp/imex-n7-composition.XXXXXX)
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
    if [ "$diagnostic" = 1 ]; then
        PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/diagnose_cases.py" \
            "$build_dir/$profile/test_imex_stages" "$profile"
    else
        PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/check_cases.py" \
            "$build_dir/$profile/test_imex_stages" "$profile" $equilibrium_only
    fi
done
if [ "$diagnostic" = 1 ]; then
    echo 'DIAGNOSTIC ONLY: N7-B remains open; inspect both diagnosis.json files.'
fi
