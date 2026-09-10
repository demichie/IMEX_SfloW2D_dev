#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-topography-faces.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -fcheck=all -fbacktrace -c \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/geometry_2d.f90" \
    "$test_dir/test_topography_faces.f90"

gfortran -O0 -g -fcheck=all -fbacktrace \
    parameters_2d.o geometry_2d.o test_topography_faces.o \
    -llapack -o test_topography_faces

./test_topography_faces
