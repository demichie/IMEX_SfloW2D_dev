#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/.." && pwd)

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/IMEX_SfloW2D" >&2
    exit 2
fi

executable=$1
case "$executable" in
    /*) ;;
    *)
        executable=$(CDPATH= cd -- "$(dirname -- "$executable")" && pwd)/$(basename -- "$executable")
        ;;
esac
if [ ! -x "$executable" ]; then
    echo "Solver executable is missing or not executable: $executable" >&2
    exit 2
fi

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/imex-integration-tests.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM

for test_name in TEST_2D TEST_BOX TEST_BW TEST_CF TEST_TAAL TEST_DAMBREAK; do
    cp -R "$test_dir/$test_name" "$work_dir/"
done

run_case() {
    case_dir=$1
    expected_output=$2
    fixture=${3:-}
    (
        cd "$case_dir"
        sh ./cleanFolder.sh
        case "$fixture" in
            cf)
                unzip -q "$repo_dir/EXAMPLES/EXAMPLE_CF/topo.zip"
                ;;
            taal)
                cp "$repo_dir/EXAMPLES/EXAMPLE_TAAL/topography_dem.asc" .
                ;;
        esac
        if ! "$executable" > solver.log 2>&1; then
            echo "FAIL: solver exited unsuccessfully in $case_dir" >&2
            cat solver.log >&2
            exit 1
        fi
        if [ ! -f "$expected_output" ]; then
            echo "FAIL: expected output was not created in $case_dir: $expected_output" >&2
            cat solver.log >&2
            exit 1
        fi
        sh ./cleanFolder.sh
    )
    echo "PASS: $(basename "$case_dir")"
}

run_case "$work_dir/TEST_2D" "TEST_2D_0001.q_2d"
run_case "$work_dir/TEST_BOX" "TEST_BOX_400_0001.q_2d"
run_case "$work_dir/TEST_BW" "TEST_BW_400_0001.q_2d"

run_case "$work_dir/TEST_CF" "TEST_CF_0001.q_2d" cf

run_case "$work_dir/TEST_TAAL" "TEST_TAAL_0001.q_2d" taal

for case_dir in "$work_dir"/TEST_DAMBREAK/*/; do
    [ -d "$case_dir" ] || continue
    run_case "$case_dir" "dambreak2D_0001.q_2d"
done

echo "PASS: legacy solver integration tests"
