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
    *)
        executable_dir=$(CDPATH= cd -- "$(dirname -- "$executable")" && pwd)
        executable="$executable_dir/$(basename -- "$executable")"
        ;;
esac

template="$repo_dir/EXAMPLES/EXAMPLE_INCLINED_EXCAVATION_2D/IMEX_SfloW2D.template"
reference="$test_dir/reference_metrics.csv"
work_dir=$(mktemp -d /tmp/imex-pccu-gate-h.XXXXXX)
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM

for mode in g1 slope_curvature
do
    for item in 0.4:dx0p4 0.2:dx0p2 0.1:dx0p1
    do
        dx=${item%%:*}
        tag=${item#*:}
        case_dir="$work_dir/$mode-$tag"
        python3 "$test_dir/generate_case.py" "$dx" "$mode" "$case_dir" "$template"
        (
            cd "$case_dir"
            OMP_NUM_THREADS=1 "$executable" > run.log
        )
        python3 "$test_dir/analyze_case.py" "$case_dir" "$mode" "$reference"
    done
done

echo "PASS: Gate-H G=1 and slope-plus-curvature excavation refinements"
echo "      match the dry-safe Python reference fingerprints"
