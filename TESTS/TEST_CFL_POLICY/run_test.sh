#!/bin/sh
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
if [ "$#" -ne 1 ]; then
    echo "Usage: sh $0 /path/to/IMEX_SfloW2D" >&2
    exit 2
fi
work_dir=$(mktemp -d /tmp/imex-cfl-policy.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/test_checker.py"
python3 "$repo_dir/TESTS/TEST_PCCU_GATE_H/build_dt_observer.py" \
    "$1" "$work_dir/observer-build" "$repo_dir" --include-spatial
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/check_cases.py" \
    "$1" "$work_dir/observer-build/observe_accepted_times" "$work_dir" "$repo_dir"
