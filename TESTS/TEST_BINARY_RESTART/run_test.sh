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

# Enable the evolving bed in both trajectories. The reference run crosses the
# same output boundary used to write the intermediate restart, so both paths
# have identical step clipping and exercise B_vertex restart fidelity.
sed -i.bak -e 's/^ T_END=.*/ T_END=  3.000000000000000E+000,/' \
    -e 's/^ TOPO_CHANGE_FLAG = F,/ TOPO_CHANGE_FLAG = T,/' \
    "$reference_dir/IMEX_sfloW2D.inp"
sed -i.bak -e 's/^ TOPO_CHANGE_FLAG = F,/ TOPO_CHANGE_FLAG = T,/' \
    "$resumed_dir/IMEX_sfloW2D.inp"
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

latest_snapshot() {
    snapshot_dir=$1
    set -- "$snapshot_dir"/TEST_2D_[0-9][0-9][0-9][0-9].q_2d
    if [ ! -f "$1" ]; then
        echo "FAIL: no conservative snapshots found in $snapshot_dir" >&2
        ls -la "$snapshot_dir" >&2
        cat "$snapshot_dir/run.log" "$snapshot_dir/resumed_run.log" 2>/dev/null || true
        exit 1
    fi
    latest_snapshot_file=
    for snapshot in "$@"; do
        latest_snapshot_file=$snapshot
    done
    printf '%s\n' "$latest_snapshot_file"
}

reference_snapshot=$(latest_snapshot "$reference_dir")
resumed_snapshot=$(latest_snapshot "$resumed_dir")
if [ "$(basename "$reference_snapshot")" != "$(basename "$resumed_snapshot")" ]; then
    echo "FAIL: full and resumed runs ended at different output indices" >&2
    echo "  reference: $reference_snapshot" >&2
    echo "  resumed:   $resumed_snapshot" >&2
    exit 1
fi
cmp "$reference_snapshot" "$resumed_snapshot"

echo "PASS: evolving-topography binary restart is bit-for-bit identical"

# Also exercise original runout references, complete/segmented CSV histories,
# read-only refusal of incompatible files and restored positive-EPS_STOP state.
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/test_runout_restart.py" "$executable"
