#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-stochastic-random.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -Wall -Wextra -fcheck=all -fbacktrace \
    "$repo_dir/src/parameters_2d.f90" \
    "$repo_dir/src/stochastic_random_2d.f90" \
    "$test_dir/test_stochastic_random.f90" \
    -o test_stochastic_random

./test_stochastic_random
