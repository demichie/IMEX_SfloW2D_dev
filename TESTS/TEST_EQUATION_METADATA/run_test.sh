#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-equation-metadata.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -Wall -fcheck=all -fbacktrace -c \
    "$repo_dir/src/equation_metadata_2d.f90" \
    "$test_dir/test_equation_metadata.f90"

gfortran -O0 -g -Wall -fcheck=all -fbacktrace \
    equation_metadata_2d.o test_equation_metadata.o \
    -o test_equation_metadata

./test_equation_metadata
