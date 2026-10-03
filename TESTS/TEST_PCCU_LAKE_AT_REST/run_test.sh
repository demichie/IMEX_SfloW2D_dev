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

work_dir=$(mktemp -d /tmp/imex-pccu-lake.XXXXXX)
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
cp "$test_dir/generate_case.py" "$work_dir/"

(
    cd "$work_dir"

    run_case() {
        label=$1
        bed_mode=$2
        shift 2
        rm -f lakeRest_0001.q_2d lakeRest_0001.1thread.q_2d
        python3 generate_case.py --bed-mode "$bed_mode" "$@"
        OMP_NUM_THREADS=1 "$executable" > "run_${label}_1_thread.log"
        cp lakeRest_0001.q_2d lakeRest_0001.1thread.q_2d
        OMP_NUM_THREADS=4 "$executable" > "run_${label}_4_threads.log"
        cmp lakeRest_0001.1thread.q_2d lakeRest_0001.q_2d
        python3 generate_case.py --check lakeRest_0001.q_2d \
            --bed-mode "$bed_mode"
        cp lakeRest_0001.q_2d "lakeRest_${label}.q_2d"
    }

    for bed_mode in planar one-cell; do
        run_case "${bed_mode}_G1" "$bed_mode"
        run_case "${bed_mode}_curvature" "$bed_mode" --curvature
        cmp "lakeRest_${bed_mode}_G1.q_2d" \
            "lakeRest_${bed_mode}_curvature.q_2d"

        run_case "${bed_mode}_slope" "$bed_mode" --slope-correction
        run_case "${bed_mode}_slope_curvature" "$bed_mode" \
            --slope-correction --curvature
        cmp "lakeRest_${bed_mode}_slope.q_2d" \
            "lakeRest_${bed_mode}_slope_curvature.q_2d"
    done
)

echo "PASS: lake at rest is thread-reproducible and unchanged for all slope/curvature combinations on planar and one-cell-ramp beds"
