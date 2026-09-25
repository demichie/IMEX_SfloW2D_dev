#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-pccu.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -fcheck=all -fbacktrace -c \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/nonconservative_2d.f90" \
    "$repo_dir/src/pccu_2d.f90" \
    "$test_dir/test_pccu.f90"

gfortran -O0 -g -fcheck=all -fbacktrace \
    parameters_2d.o nonconservative_2d.o pccu_2d.o test_pccu.o \
    -o test_pccu

./test_pccu
