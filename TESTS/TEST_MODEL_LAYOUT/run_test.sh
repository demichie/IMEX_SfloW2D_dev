#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-model-layout.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"

gfortran -O0 -g -Wall -Wextra -fcheck=all -fbacktrace -c \
    "$repo_dir/src/model_layout_2d.f90" \
    "$test_dir/test_model_layout.f90" \
    "$test_dir/test_model_layout_reject.f90"

gfortran -O0 -g -Wall -Wextra -fcheck=all -fbacktrace \
    model_layout_2d.o test_model_layout.o -o test_model_layout

./test_model_layout

gfortran -O0 -g -Wall -Wextra -fcheck=all -fbacktrace \
    model_layout_2d.o test_model_layout_reject.o \
    -o test_model_layout_reject

if ./test_model_layout_reject >reject.log 2>&1; then
    echo "FAIL: N_LAYERS=2 was accepted" >&2
    exit 1
fi

grep -q "Only N_LAYERS=1 is currently supported" reject.log
echo "PASS: unsupported N_LAYERS=2 rejected explicitly"
