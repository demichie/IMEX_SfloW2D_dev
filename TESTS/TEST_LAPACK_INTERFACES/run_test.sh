#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-lapack-interfaces.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/lapack_interfaces_2d.f90" \
    "$test_dir/test_lapack_interfaces.f90" \
    -llapack -o test_lapack_interfaces

./test_lapack_interfaces
