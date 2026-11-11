#!/bin/sh
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-layered-storage.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
cd "$build_dir"
gfortran -O0 -g -Wall -Wextra -fcheck=all -fbacktrace \
    "$repo_dir/src/core/parameters_2d.f90" \
    "$repo_dir/src/core/model_layout_2d.f90" \
    "$repo_dir/src/state/layered_fluid_storage_2d.f90" \
    "$test_dir/test_layered_fluid_storage.f90" \
    -o test_layered_fluid_storage
./test_layered_fluid_storage
