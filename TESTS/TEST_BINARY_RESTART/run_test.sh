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

prepare_solver_input() {
    case_dir=$1
    input_tmp="$case_dir/IMEX_SfloW2D.inp.new"
    cat "$case_dir/IMEX_sfloW2D.inp" > "$input_tmp"
    rm "$case_dir/IMEX_sfloW2D.inp"
    mv "$input_tmp" "$case_dir/IMEX_SfloW2D.inp"
}
prepare_solver_input "$reference_dir"
prepare_solver_input "$resumed_dir"

# Enable the evolving bed in both trajectories. The reference run crosses the
# same output boundary used to write the intermediate restart, so both paths
# have identical step clipping and exercise B_vertex restart fidelity.
sed -i.bak -e 's/^ T_END=.*/ T_END=  3.000000000000000E+000,/' \
    -e 's/^ TOPO_CHANGE_FLAG = F,/ TOPO_CHANGE_FLAG = T,/' \
    "$reference_dir/IMEX_SfloW2D.inp"
sed -i.bak -e 's/^ TOPO_CHANGE_FLAG = F,/ TOPO_CHANGE_FLAG = T,/' \
    "$resumed_dir/IMEX_SfloW2D.inp"
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
        IMEX_SfloW2D.inp
    "$executable" > resumed_run.log
)

cmp "$reference_dir/OUTPUT/TEST_2D/restart/restart.bin" \
    "$resumed_dir/OUTPUT/TEST_2D/restart/restart.bin"

latest_snapshot() {
    snapshot_dir=$1
    set -- "$snapshot_dir"/OUTPUT/TEST_2D/snapshots/*/TEST_2D_[0-9][0-9][0-9][0-9].q_2d
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

# An external checkpoint may use a long absolute path containing spaces. It
# starts a fresh output run but must not authorize reuse of an existing run.
external_case="$work_dir/external-restart-case"
external_checkpoint_dir="$work_dir/external restart checkpoint directory with a path longer than forty characters"
external_checkpoint="$external_checkpoint_dir/saved restart checkpoint.bin"
external_output_base="external output base with spaces"
mkdir -p "$external_checkpoint_dir" "$external_case"
cp "$reference_dir/OUTPUT/TEST_2D/restart/restart.bin" "$external_checkpoint"
cp -R "$repo_dir/TESTS/TEST_2D/." "$external_case/"
prepare_solver_input "$external_case"
[ "${#external_checkpoint}" -gt 40 ] || {
    echo "FAIL: external checkpoint fixture path is not longer than 40 characters" >&2
    exit 1
}
(
    cd "$external_case"
    awk -v checkpoint="$external_checkpoint" -v output_base="$external_output_base" '
        /RESTART_FILES=/ {
            print " RESTART_FILES=\"" checkpoint "\","
            next
        }
        /^ T_END=/ {
            print " T_END=  4.000000000000000E+000,"
            next
        }
        /RUN_NAME=/ {
            print
            printf " OUTPUT_BASE_DIRECTORY=\"%s\",\n", output_base
            next
        }
        /^ TOPO_CHANGE_FLAG =/ {
            print " TOPO_CHANGE_FLAG = T,"
            next
        }
        { print }
    ' IMEX_SfloW2D.inp > external.inp
    mv external.inp IMEX_SfloW2D.inp
    "$executable" > external-run.log 2>&1 || {
        echo "FAIL: external restart with a long path containing spaces failed" >&2
        cat external-run.log >&2
        exit 1
    }
    run_output="$external_output_base/TEST_2D"
    [ -f "$run_output/restart/restart.bin" ] || {
        echo "FAIL: external restart did not write into OUTPUT_BASE_DIRECTORY" >&2
        exit 1
    }
    manifest_checksum=$(cksum "$run_output/metadata/output_manifest.csv")
    restart_checksum=$(cksum "$run_output/restart/restart.bin")

    if "$executable" > external-collision.log 2>&1; then
        echo "FAIL: external checkpoint implicitly authorized existing output reuse" >&2
        exit 1
    fi
    grep -q 'run output directory already exists' external-collision.log || {
        echo "FAIL: external restart collision did not report the expected error" >&2
        cat external-collision.log >&2
        exit 1
    }
    [ "$manifest_checksum" = "$(cksum "$run_output/metadata/output_manifest.csv")" ] || {
        echo "FAIL: rejected external restart changed the existing manifest" >&2
        exit 1
    }
    [ "$restart_checksum" = "$(cksum "$run_output/restart/restart.bin")" ] || {
        echo "FAIL: rejected external restart changed the existing checkpoint" >&2
        exit 1
    }

    awk '{ print; if ($0 ~ /RUN_NAME/) print "  OUTPUT_OVERWRITE=T," }' \
        IMEX_SfloW2D.inp > overwrite.inp
    mv overwrite.inp IMEX_SfloW2D.inp
    output_checksums=$(find "$run_output" -type f -exec cksum {} \; | LC_ALL=C sort)
    if "$executable" > external-stale-overwrite.log 2>&1; then
        echo "FAIL: OUTPUT_OVERWRITE allowed an external checkpoint to duplicate manifest output" >&2
        cat external-stale-overwrite.log >&2
        exit 1
    fi
    grep -q 'manifest' external-stale-overwrite.log || {
        echo "FAIL: stale external restart did not report the manifest conflict" >&2
        cat external-stale-overwrite.log >&2
        exit 1
    }
    updated_output_checksums=$(find "$run_output" -type f -exec cksum {} \; | LC_ALL=C sort)
    [ "$output_checksums" = "$updated_output_checksums" ] || {
        echo "FAIL: rejected stale external restart changed existing output files" >&2
        printf '%s\n' "$output_checksums" > output-checksums-before.txt
        printf '%s\n' "$updated_output_checksums" > output-checksums-after.txt
        diff -u output-checksums-before.txt output-checksums-after.txt >&2 || true
        exit 1
    }

    continuation_checkpoint="$external_checkpoint_dir/continuation restart checkpoint.bin"
    cp "$run_output/restart/restart.bin" "$continuation_checkpoint"
    awk -v checkpoint="$continuation_checkpoint" '
        /RESTART_FILES=/ {
            print " RESTART_FILES=\"" checkpoint "\","
            next
        }
        /^ T_END=/ {
            print " T_END=  5.500000000000000E+000,"
            next
        }
        { print }
    ' IMEX_SfloW2D.inp > continuation.inp
    mv continuation.inp IMEX_SfloW2D.inp
    "$executable" > external-overwrite.log 2>&1 || {
        echo "FAIL: explicit overwrite did not allow compatible external restart continuation" >&2
        cat external-overwrite.log >&2
        exit 1
    }
)
echo "PASS: long external restart paths and explicit output reuse policy"

# Also exercise original runout references, complete/segmented CSV histories,
# read-only refusal of incompatible files and restored positive-EPS_STOP state.
PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/test_runout_restart.py" "$executable"
