#!/usr/bin/env bash
# Exercise the compiled solver in the exact image, not just its dynamic links.
set -euo pipefail
if (( $# != 1 )); then
    echo "Usage: $0 IMAGE_REFERENCE" >&2
    exit 2
fi
image=$1
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

cp TESTS/TEST_2D/IMEX_sfloW2D.inp "$work/IMEX_SfloW2D.inp"
cp TESTS/TEST_2D/INIT_2D_0000.q_2d "$work/"
cp TESTS/TEST_2D/topography_dem.asc "$work/"

docker run --rm --entrypoint /bin/sh "$image" -ec \
    'test -x /usr/local/bin/IMEX_SfloW2D && ! ldd /usr/local/bin/IMEX_SfloW2D | grep -q "not found"'

if ! timeout --signal=TERM --kill-after=10s 600s docker run --rm \
    --user "$(id -u):$(id -g)" \
    --mount "type=bind,src=$work,dst=/work" --workdir /work \
    "$image" > "$work/solver.log" 2>&1; then
    tail -n 100 "$work/solver.log" >&2
    echo "FAIL: container did not finish the TEST_2D solver fixture" >&2
    exit 1
fi

manifest="$work/OUTPUT/TEST_2D/metadata/output_manifest.csv"
test -s "$manifest" || {
    tail -n 100 "$work/solver.log" >&2
    echo "FAIL: container produced no output manifest" >&2
    exit 1
}
grep -q '^snapshot_index,time_seconds$' "$manifest"
if ! find "$work/OUTPUT/TEST_2D/snapshots" -type f \
    -name '*.q_2d' -print -quit | grep -q .; then
    echo "FAIL: container produced no conservative-state snapshots" >&2
    exit 1
fi
echo "PASS: container linkage and TEST_2D end-to-end solver output"
