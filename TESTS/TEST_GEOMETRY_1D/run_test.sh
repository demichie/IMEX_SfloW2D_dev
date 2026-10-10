#!/bin/sh
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-geometry-1d.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $build_dir"
else
    trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
fi
for profile in strict optimized; do
    mkdir "$build_dir/$profile"
    cd "$build_dir/$profile"
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all' ;;
        optimized) flags='-Ofast -funroll-all-loops' ;;
    esac
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -c \
        "$repo_dir/src/parameters_2d.f90" "$repo_dir/src/geometry_2d.f90" \
        "$test_dir/test_geometry_1d.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        parameters_2d.o geometry_2d.o test_geometry_1d.o -llapack -o test_geometry_1d
    PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/check_cases.py" "$build_dir/$profile/test_geometry_1d" "$profile" "$@"
done
