#!/bin/sh
# Source mass balances and invalid reduced-gravity input rejection.
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)

if [ "$#" -ge 1 ]; then
    executable=$1
elif [ -x "$repo_dir/src/IMEX_SfloW2D" ]; then
    executable="$repo_dir/src/IMEX_SfloW2D"
else
    echo "Usage: $0 /path/to/IMEX_SfloW2D [tolerance]" >&2
    exit 2
fi
tolerance=${2:-1.0e-9}

case "$executable" in
    /*) ;;
    *)
        executable="$(CDPATH= cd -- "$(dirname -- "$executable")" && pwd)/$(basename -- "$executable")"
        ;;
esac

work_dir=$(mktemp -d /tmp/imex-source-mass.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

python3 -B "$test_dir/make_cases.py" "$work_dir"

failed=0
for case_dir in "$work_dir"/*/; do
    case_name=$(basename "$case_dir")
    if [ "$case_name" = source_negative_reduced_gravity ]; then
        if ( cd "$case_dir" && "$executable" > run.log 2>&1 ); then
            echo "FAIL: solver accepted a source with non-positive reduced gravity"
            failed=1
            continue
        fi
        if ! grep -q 'Source reduced gravity must be positive when RI_SOURCE is specified' \
                "$case_dir/run.log"; then
            echo "FAIL: invalid source did not report the expected validation error"
            tail -20 "$case_dir/run.log"
            failed=1
            continue
        fi
        if [ -e "$case_dir/OUTPUT" ]; then
            echo "FAIL: invalid source created solver output before rejection"
            failed=1
            continue
        fi
        echo "PASS: rejected source with non-positive reduced gravity"
        continue
    fi

    if ! ( cd "$case_dir" && IMEX_SOURCE_BOUNDARY_BUDGET=1 "$executable" > run.log 2>&1 ); then
        echo "FAIL: solver exited unsuccessfully in $case_dir"
        tail -20 "$case_dir/run.log"
        failed=1
        continue
    fi
    if ! python3 -B "$test_dir/check_mass.py" "$case_dir" "$tolerance"; then
        case "$case_name" in
            lateral_boundary_*)
                echo "Source initialization diagnostics for $case_name:"
                grep -E "Source (density|reduced gravity|Richardson number|velocity|thickness|mass flow rate)" "$case_dir/run.log" || :
                echo "Boundary flux budget near the final time:"
                grep "SOURCE_BOUNDARY_BUDGET" "$case_dir/run.log" | tail -8 || :
                echo "Final solver diagnostic lines:"
                tail -12 "$case_dir/run.log"
                ;;
        esac
        failed=1
    fi
done

if [ "$failed" -ne 0 ]; then
    echo "FAIL: source mass injection"
    exit 1
fi
echo "PASS: source mass injection and validation"
