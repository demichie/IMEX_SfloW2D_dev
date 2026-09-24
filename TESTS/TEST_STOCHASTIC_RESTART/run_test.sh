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

work_dir=$(mktemp -d /tmp/imex-stochastic-restart.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

prepare_case() {
    case_dir=$1
    mkdir "$case_dir"
    cp -R "$repo_dir/TESTS/TEST_2D/." "$case_dir/"

    # The stochastic transport equation adds one conservative variable. Append
    # its initial hZ value to each non-empty record of the text initial state.
    awk 'NF { print $0, "0.000000000000E+00"; next } { print }' \
        "$case_dir/INIT_2D_0000.q_2d" > \
        "$case_dir/INIT_STOCHASTIC_0000.q_2d"

    awk '
        /^ RESTART_FILES=/ {
            print " RESTART_FILES=\"INIT_STOCHASTIC_0000.q_2d\","; next
        }
        /^ LIQUID_FLAG = F,/ {
            print
            print " STOCHASTIC_FLAG = T,"
            next
        }
        /^ LIMITER=/ {
            print " LIMITER= 1 1 1 1 1 1,"
            next
        }
        /^ T_BCW%VALUE=/ {
            print
            print " STOCH_BCW%FLAG=1,"
            print " STOCH_BCW%VALUE=0.0D0,"
            next
        }
        /^ T_BCE%VALUE=/ {
            print
            print " STOCH_BCE%FLAG=1,"
            print " STOCH_BCE%VALUE=0.0D0,"
            next
        }
        /^ T_BCS%VALUE=/ {
            print
            print " STOCH_BCS%FLAG=1,"
            print " STOCH_BCS%VALUE=0.0D0,"
            next
        }
        /^ T_BCN%VALUE=/ {
            print
            print " STOCH_BCN%FLAG=1,"
            print " STOCH_BCN%VALUE=0.0D0,"
            next
        }
        /^ RHEOLOGY_MODEL =/ {
            print " RHEOLOGY_MODEL = 9,"
            print " MU_0 = 0.10D0,"
            print " MU_INF = 0.50D0,"
            print " FR_0 = 1.0D0,"
            next
        }
        /^ FRICTION_FACTOR =/ { next }
        /^&SOLID_TRANSPORT_PARAMETERS/ {
            print "&STOCHASTIC_PARAMETERS"
            print " SYM_NOISE = 0.0D0,"
            print " STD_MAX = 0.05D0,"
            print " TAU_STOCHASTIC = 1.0D0,"
            print " LENGTH_SPATIAL_CORR = 0.0D0,"
            print " STOCH_TRANSPORT_FLAG = T,"
            print " STOCHASTIC_SEED = 24680,"
            print " /"
            print ""
            print $0
            next
        }
        { print }
    ' "$case_dir/IMEX_sfloW2D.inp" > "$case_dir/IMEX_sfloW2D.inp.new"
    mv "$case_dir/IMEX_sfloW2D.inp.new" "$case_dir/IMEX_sfloW2D.inp"
}

reference_dir="$work_dir/reference"
resumed_dir="$work_dir/resumed"
prepare_case "$reference_dir"
prepare_case "$resumed_dir"

# Both trajectories cross the same output boundary so timestep clipping is
# identical. The second trajectory is interrupted and resumed at that boundary.
sed -i.bak -e 's/^ T_END=.*/ T_END=  3.000000000000000E+000,/' \
    "$reference_dir/IMEX_sfloW2D.inp"
(
    cd "$reference_dir"
    "$executable" > run.log
)

(
    cd "$resumed_dir"
    "$executable" > first_run.log
    sed -i.bak \
        -e 's/^ T_END=.*/ T_END=  3.000000000000000E+000,/' \
        -e 's/^ RESTART_FILES=.*/ RESTART_FILES="restart.bin",/' \
        IMEX_sfloW2D.inp
    "$executable" > resumed_run.log
)

if [ ! -f "$reference_dir/restart.bin" ] || \
   [ ! -f "$resumed_dir/restart.bin" ]; then
    echo "FAIL: stochastic run did not produce the expected restart file" >&2
    tail -80 "$reference_dir/run.log" >&2
    tail -80 "$resumed_dir/first_run.log" >&2
    exit 1
fi

cmp "$reference_dir/restart.bin" "$resumed_dir/restart.bin"
cmp "$reference_dir/TEST_2D_0002.q_2d" \
    "$resumed_dir/TEST_2D_0002.q_2d"

# Columns 1-2 are coordinates and columns 3-8 are the six conservative
# variables. Confirm that the OU/transport variable was actually evolved.
if ! awk 'NF && $8 != 0.0 { found=1 } END { exit !found }' \
    "$reference_dir/TEST_2D_0002.q_2d"; then
    echo "FAIL: the stochastic conservative variable remained identically zero" >&2
    exit 1
fi

echo "PASS: stochastic transport and restarted OU evolution are bit-for-bit identical"
