#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/.." && pwd)

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 /path/to/IMEX_SfloW2D [TEST_DYNAMIC_FINH]" >&2
    exit 2
fi

executable=$1
selected_case=${2:-}
case "$selected_case" in
    ""|TEST_DYNAMIC_FINH) ;;
    *)
        echo "Unsupported focused integration case: $selected_case" >&2
        exit 2
        ;;
esac
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

prepare_solver_input() {
    fixture_dir=$1
    if [ -f "$fixture_dir/IMEX_sfloW2D.inp" ]; then
        mv "$fixture_dir/IMEX_sfloW2D.inp" "$fixture_dir/IMEX_input.tmp"
        mv "$fixture_dir/IMEX_input.tmp" "$fixture_dir/IMEX_SfloW2D.inp"
    fi
}

copy_writable_fixture() {
    cp -R "$1" "$2"
    chmod -R u+w "$2"
    prepare_solver_input "$2"
}

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/imex-integration-tests.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM

if [ -n "$selected_case" ]; then
    copy_writable_fixture "$test_dir/TEST_DAMBREAK/$selected_case" \
        "$work_dir/$selected_case"
else
    for test_name in TEST_2D TEST_BOX TEST_BW TEST_CF TEST_TAAL TEST_DAMBREAK; do
        copy_writable_fixture "$test_dir/$test_name" "$work_dir/$test_name"
    done
fi

run_case() {
    case_dir=$1
    run_name=$2
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
        solver_log=solver.log
        if [ -n "${IMEX_TEST_LOG_DIR:-}" ]; then
            mkdir -p "$IMEX_TEST_LOG_DIR"
            case_name=$(basename "$case_dir")
            solver_log="$IMEX_TEST_LOG_DIR/$case_name.solver.log"
        fi

        solver_status=0
        if [ -n "${IMEX_TEST_LOG_DIR:-}" ]; then
            run_started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
            run_started_epoch=$(date +%s)
        fi
        if [ -n "${IMEX_TEST_CASE_TIMEOUT:-}" ]; then
            if timeout --signal=TERM --kill-after=30s \
                "$IMEX_TEST_CASE_TIMEOUT" "$executable" > "$solver_log" 2>&1; then
                :
            else
                solver_status=$?
            fi
        elif "$executable" > "$solver_log" 2>&1; then
            :
        else
            solver_status=$?
        fi
        if [ -n "${IMEX_TEST_LOG_DIR:-}" ]; then
            run_finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
            run_finished_epoch=$(date +%s)
            printf 'started_utc=%s\nfinished_utc=%s\nelapsed_seconds=%s\n' \
                "$run_started_utc" "$run_finished_utc" \
                "$((run_finished_epoch - run_started_epoch))" \
                > "$IMEX_TEST_LOG_DIR/$case_name.timing.txt"
        fi
        if [ "$solver_status" -ne 0 ]; then
            echo "FAIL: solver exited unsuccessfully in $case_dir" >&2
            cat "$solver_log" >&2
            exit "$solver_status"
        fi
        set -- "$case_dir/OUTPUT/$run_name"/snapshots/*/"$run_name"_*.q_2d
        if [ ! -f "$1" ]; then
            echo "FAIL: expected conservative snapshot was not created for $run_name" >&2
            find "$case_dir/OUTPUT" -maxdepth 5 -type f -print >&2 || true
            cat "$solver_log" >&2
            exit 1
        fi
        [ -f "$case_dir/OUTPUT/$run_name/metadata/$run_name.bak" ] || {
            echo "FAIL: input backup is missing from per-run metadata for $run_name" >&2
            exit 1
        }
        grep -q '^snapshot_index,time_seconds$' \
            "$case_dir/OUTPUT/$run_name/metadata/output_manifest.csv" || {
            echo "FAIL: snapshot manifest is missing or malformed for $run_name" >&2
            exit 1
        }
        grep -Eq '^[0-9]+, *[-+]?[0-9]' \
            "$case_dir/OUTPUT/$run_name/metadata/output_manifest.csv" || {
            echo "FAIL: snapshot manifest contains no index/time records for $run_name" >&2
            exit 1
        }
        if find "$case_dir" -maxdepth 1 -type f \
            \( -name "$run_name*.q_2d" -o -name "$run_name*.nc" -o \
               -name "$run_name*.csv" -o -name "$run_name.bak" -o \
               -name "restart.bin" -o -name "dem_esri*.asc" \) \
            | grep -q .; then
            echo "FAIL: solver output escaped the per-run directory for $run_name" >&2
            exit 1
        fi
        if [ "$run_name" = TEST_2D ]; then
            before=$(cksum "$1")
            if "$executable" > collision.log 2>&1; then
                echo "FAIL: an existing run directory was accepted without OUTPUT_OVERWRITE" >&2
                exit 1
            fi
            grep -q 'run output directory already exists' collision.log || {
                echo "FAIL: output collision did not report the expected error" >&2
                cat collision.log >&2
                exit 1
            }
            [ "$before" = "$(cksum "$1")" ] || {
                echo "FAIL: rejected output collision changed an existing snapshot" >&2
                exit 1
            }
            awk '{ print; if ($0 ~ /RUN_NAME/) print "  OUTPUT_OVERWRITE=T," }' \
                IMEX_SfloW2D.inp > IMEX_SfloW2D.inp.tmp
            mv IMEX_SfloW2D.inp.tmp IMEX_SfloW2D.inp
            "$executable" > overwrite.log 2>&1 || {
                echo "FAIL: OUTPUT_OVERWRITE=T did not allow an existing run directory" >&2
                cat overwrite.log >&2
                exit 1
            }
        fi
        sh ./cleanFolder.sh
    )
    echo "PASS: $(basename "$case_dir")"
}

if [ -n "$selected_case" ]; then
    run_case "$work_dir/$selected_case" "dambreak2D"
    exit 0
fi

check_output_path_safety() {
    safety_dir=$work_dir/output-path-safety
    copy_writable_fixture "$test_dir/TEST_2D" "$safety_dir"
    (
        cd "$safety_dir"
        sh ./cleanFolder.sh
        sed 's/RUN_NAME="TEST_2D"/RUN_NAME="..\/escape"/' \
            IMEX_SfloW2D.inp > invalid.inp
        mv invalid.inp IMEX_SfloW2D.inp
        if "$executable" > invalid-name.log 2>&1; then
            echo "FAIL: path traversal in RUN_NAME was accepted" >&2
            exit 1
        fi
        grep -q 'RUN_NAME may contain only' invalid-name.log || {
            echo "FAIL: invalid RUN_NAME did not report the expected error" >&2
            cat invalid-name.log >&2
            exit 1
        }
        [ ! -e OUTPUT ] || {
            echo "FAIL: invalid RUN_NAME created an output directory" >&2
            exit 1
        }
    )

    copy_writable_fixture "$test_dir/TEST_2D" "$safety_dir/invalid-name-existing-case"
    (
        cd "$safety_dir/invalid-name-existing-case"
        sh ./cleanFolder.sh
        "$executable" > valid-name.log 2>&1 || {
            echo "FAIL: could not prepare existing output for invalid-input check" >&2
            cat valid-name.log >&2
            exit 1
        }
        before_files=$(find OUTPUT -type f -exec cksum {} \; | LC_ALL=C sort)
        before_dirs=$(find OUTPUT -type d -print | LC_ALL=C sort)
        sed 's/RUN_NAME="TEST_2D"/RUN_NAME="..\/escape"/' \
            IMEX_SfloW2D.inp > invalid.inp
        mv invalid.inp IMEX_SfloW2D.inp
        if "$executable" > invalid-name-existing.log 2>&1; then
            echo "FAIL: path traversal in RUN_NAME was accepted with existing output" >&2
            exit 1
        fi
        grep -q 'RUN_NAME may contain only' invalid-name-existing.log || {
            echo "FAIL: invalid RUN_NAME did not report the expected error" >&2
            cat invalid-name-existing.log >&2
            exit 1
        }
        after_files=$(find OUTPUT -type f -exec cksum {} \; | LC_ALL=C sort)
        after_dirs=$(find OUTPUT -type d -print | LC_ALL=C sort)
        [ "$before_files" = "$after_files" ] && [ "$before_dirs" = "$after_dirs" ] || {
            echo "FAIL: invalid RUN_NAME changed the existing output tree" >&2
            exit 1
        }
        [ ! -e escape ] || {
            echo "FAIL: invalid RUN_NAME created output outside OUTPUT" >&2
            exit 1
        }
    )

    copy_writable_fixture "$test_dir/TEST_2D" "$safety_dir/symlink-case"
    (
        cd "$safety_dir/symlink-case"
        sh ./cleanFolder.sh
        mkdir symlink-target
        ln -s symlink-target unsafe-base
        awk '{ print; if ($0 ~ /RUN_NAME/) print "  OUTPUT_BASE_DIRECTORY=\"unsafe-base\"," }' \
            IMEX_SfloW2D.inp > symlink.inp
        mv symlink.inp IMEX_SfloW2D.inp
        if "$executable" > symlink.log 2>&1; then
            echo "FAIL: symlink in OUTPUT_BASE_DIRECTORY was accepted" >&2
            exit 1
        fi
        [ -z "$(find symlink-target -mindepth 1 -print -quit)" ] || {
            echo "FAIL: symlink target was modified during path validation" >&2
            exit 1
        }
    )

    copy_writable_fixture "$test_dir/TEST_2D" "$safety_dir/symlink-file-case"
    (
        cd "$safety_dir/symlink-file-case"
        sh ./cleanFolder.sh
        mkdir -p OUTPUT/TEST_2D/metadata
        printf 'protected sentinel\n' > "$safety_dir/sentinel.txt"
        ln -s "$safety_dir/sentinel.txt" OUTPUT/TEST_2D/metadata/TEST_2D.bak
        awk '{ print; if ($0 ~ /RUN_NAME/) print "  OUTPUT_OVERWRITE=T," }' \
            IMEX_SfloW2D.inp > overwrite.inp
        mv overwrite.inp IMEX_SfloW2D.inp
        if "$executable" > symlink-file.log 2>&1; then
            echo "FAIL: symlink output file was followed during overwrite" >&2
            exit 1
        fi
        grep -q 'existing run output contains a symlink' symlink-file.log || {
            echo "FAIL: symlink output file did not report the expected error" >&2
            cat symlink-file.log >&2
            exit 1
        }
        [ "$(cat "$safety_dir/sentinel.txt")" = "protected sentinel" ] || {
            echo "FAIL: output symlink modified the external sentinel" >&2
            exit 1
        }
        rm OUTPUT/TEST_2D/metadata/TEST_2D.bak
        ln "$safety_dir/sentinel.txt" OUTPUT/TEST_2D/metadata/TEST_2D.bak
        if "$executable" > hardlink-file.log 2>&1; then
            echo "FAIL: hard-linked output file was accepted during overwrite" >&2
            exit 1
        fi
        tr -d '[:space:]' < hardlink-file.log | grep -q 'multiply-linkedfile' || {
            echo "FAIL: hard-linked output file did not report the expected error" >&2
            cat hardlink-file.log >&2
            exit 1
        }
        [ "$(cat "$safety_dir/sentinel.txt")" = "protected sentinel" ] || {
            echo "FAIL: output hard link modified the external sentinel" >&2
            exit 1
        }
    )

    copy_writable_fixture "$test_dir/TEST_2D" "$safety_dir/traversal-base-case"
    (
        cd "$safety_dir/traversal-base-case"
        sh ./cleanFolder.sh
        awk '{ print; if ($0 ~ /RUN_NAME/) print "  OUTPUT_BASE_DIRECTORY=\"../outside\"," }' \
            IMEX_SfloW2D.inp > traversal-base.inp
        mv traversal-base.inp IMEX_SfloW2D.inp
        if "$executable" > traversal-base.log 2>&1; then
            echo "FAIL: parent-directory component in OUTPUT_BASE_DIRECTORY was accepted" >&2
            exit 1
        fi
        [ ! -e "$work_dir/outside" ] || {
            echo "FAIL: parent-directory output path created a directory outside the test case" >&2
            exit 1
        }
    )

    copy_writable_fixture "$test_dir/TEST_2D" "$safety_dir/custom-base-case"
    (
        cd "$safety_dir/custom-base-case"
        sh ./cleanFolder.sh
        awk '{ print; if ($0 ~ /RUN_NAME/) print "  OUTPUT_BASE_DIRECTORY=\"custom output\"," }' \
            IMEX_SfloW2D.inp > custom-base.inp
        mv custom-base.inp IMEX_SfloW2D.inp
        "$executable" > custom-base.log 2>&1 || {
            echo "FAIL: valid custom output base was rejected" >&2
            cat custom-base.log >&2
            exit 1
        }
        [ -f "custom output/TEST_2D/metadata/output_manifest.csv" ] || {
            echo "FAIL: solver ignored OUTPUT_BASE_DIRECTORY" >&2
            exit 1
        }
    )
    copy_writable_fixture "$test_dir/TEST_2D" "$safety_dir/absolute-base-case"
    physical_safety_dir=$(CDPATH= cd -- "$safety_dir" && pwd -P)
    (
        cd "$safety_dir/absolute-base-case"
        sh ./cleanFolder.sh
        awk -v output="$physical_safety_dir/absolute output" \
            '{ print; if ($0 ~ /RUN_NAME/) print "  OUTPUT_BASE_DIRECTORY=\"" output "\"," }' \
            IMEX_SfloW2D.inp > absolute-base.inp
        mv absolute-base.inp IMEX_SfloW2D.inp
        "$executable" > absolute-base.log 2>&1 || {
            echo "FAIL: valid absolute OUTPUT_BASE_DIRECTORY was rejected" >&2
            cat absolute-base.log >&2
            exit 1
        }
        [ -f "$physical_safety_dir/absolute output/TEST_2D/metadata/output_manifest.csv" ] || {
            echo "FAIL: solver ignored absolute OUTPUT_BASE_DIRECTORY" >&2
            exit 1
        }
    )
    echo "PASS: output path safety"
}

run_case "$work_dir/TEST_2D" "TEST_2D"
run_case "$work_dir/TEST_BOX" "TEST_BOX_400"
run_case "$work_dir/TEST_BW" "TEST_BW_400"

run_case "$work_dir/TEST_CF" "TEST_CF" cf

run_case "$work_dir/TEST_TAAL" "TEST_TAAL" taal
check_output_path_safety

for case_dir in "$work_dir"/TEST_DAMBREAK/*/; do
    [ -d "$case_dir" ] || continue
    run_case "$case_dir" "dambreak2D"
done

echo "PASS: legacy solver integration tests"
