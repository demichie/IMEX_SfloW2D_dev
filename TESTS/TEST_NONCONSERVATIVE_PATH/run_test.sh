#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-nonconservative-path.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -Wall -Wno-unused-dummy-argument -fcheck=all -fbacktrace -c \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/nonconservative_2d.f90" \
    "$test_dir/test_nonconservative_path.f90"

gfortran -O0 -g -Wall -Wno-unused-dummy-argument -fcheck=all -fbacktrace \
    parameters_2d.o nonconservative_2d.o test_nonconservative_path.o \
    -o test_nonconservative_path

./test_nonconservative_path
