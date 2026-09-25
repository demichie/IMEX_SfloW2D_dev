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
    python3 generate_case.py
    "$executable" > run.log
    python3 generate_case.py --check lakeRest_0001.q_2d
)

echo "PASS: end-to-end PCCU lake at rest preserved to roundoff"
