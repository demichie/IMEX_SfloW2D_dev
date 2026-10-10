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
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/test_checker.py"
OMP_NUM_THREADS=${IMEX_AUDIT_FORCE_THREADS:-${OMP_NUM_THREADS:-1}}
OMP_DYNAMIC=FALSE
export OMP_NUM_THREADS OMP_DYNAMIC
python3 "$test_dir/build_dt_observer.py" "$executable" "$work_dir/observer-build" "$repo_dir"
observer="$work_dir/observer-build/observe_accepted_times"

failed=0
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
            "$executable" > run.log
        )
        python3 "$test_dir/analyze_case.py" "$case_dir" "$mode" "$reference"
        python3 "$test_dir/capture_times.py" "$observer" "$case_dir"
        if ! PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/compare_fields.py" "$case_dir" "$mode"; then
            failed=1
        fi
    done
done

if [ "$failed" -ne 0 ]; then
    echo "FAIL: Gate-H full-field evidence retained for every completed grid/mode" >&2
    exit "$failed"
fi

echo "PASS: Gate-H G=1 and slope-plus-curvature excavation refinements"
echo "      match the dry-safe Python reference fingerprints"
