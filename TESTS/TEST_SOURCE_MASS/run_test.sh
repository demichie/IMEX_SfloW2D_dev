#!/bin/sh
# Mass injected by every source type against the exact value.
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
    if ! ( cd "$case_dir" && "$executable" > run.log 2>&1 ); then
        echo "FAIL: solver exited unsuccessfully in $case_dir"
        tail -20 "$case_dir/run.log"
        failed=1
        continue
    fi
    python3 -B "$test_dir/check_mass.py" "$case_dir" "$tolerance" || failed=1
done

if [ "$failed" -ne 0 ]; then
    echo "FAIL: source mass injection"
    exit 1
fi
echo "PASS: source mass injection"
