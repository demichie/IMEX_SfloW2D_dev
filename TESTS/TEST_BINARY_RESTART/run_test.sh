#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)

if [ "$#" -eq 1 ]; then
    executable=$1
elif [ -x "$repo_dir/src/IMEX_SfloW2D" ]; then
    executable="$repo_dir/src/IMEX_SfloW2D"
elif [ -x "$repo_dir/bin/IMEX_SfloW2D" ]; then
    executable="$repo_dir/bin/IMEX_SfloW2D"
else
    echo "Usage: $0 /path/to/IMEX_SfloW2D" >&2
    exit 2
fi

case "$executable" in
    /*) ;;
    *) executable=$(CDPATH= cd -- "$(dirname -- "$executable")" && pwd)/$(basename -- "$executable") ;;
esac

work_dir=$(mktemp -d /tmp/imex-binary-restart.XXXXXX)
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM

reference_dir="$work_dir/reference"
resumed_dir="$work_dir/resumed"
mkdir "$reference_dir" "$resumed_dir"
cp -R "$repo_dir/TESTS/TEST_2D/." "$reference_dir/"
cp -R "$repo_dir/TESTS/TEST_2D/." "$resumed_dir/"

# The reference run crosses the same output boundary used to write the
# intermediate restart, so both trajectories have identical step clipping.
sed -i.bak -e 's/^ T_END=.*/ T_END=  3.000000000000000E+000,/' \
    "$reference_dir/IMEX_sfloW2D.inp"
(
    cd "$reference_dir"
    "$executable" > run.log
)

# Stop at the first output, then continue from the binary restart to t=3.
(
    cd "$resumed_dir"
    "$executable" > first_run.log
    sed -i.bak \
        -e 's/^ T_END=.*/ T_END=  3.000000000000000E+000,/' \
        -e 's/^ RESTART_FILES=.*/ RESTART_FILES="restart.bin",/' \
        IMEX_sfloW2D.inp
    "$executable" > resumed_run.log
)

cmp "$reference_dir/restart.bin" "$resumed_dir/restart.bin"
cmp "$reference_dir/TEST_2D_0002.q_2d" \
    "$resumed_dir/TEST_2D_0002.q_2d"

echo "PASS: uninterrupted and binary-restarted runs are bit-for-bit identical"
