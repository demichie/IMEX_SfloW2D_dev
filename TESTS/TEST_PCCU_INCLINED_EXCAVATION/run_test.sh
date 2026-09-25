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

example_dir="$repo_dir/EXAMPLES/EXAMPLE_INCLINED_EXCAVATION_2D"
work_dir=$(mktemp -d /tmp/imex-pccu-excavation.XXXXXX)
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
cp "$example_dir/create_example.py" "$example_dir/IMEX_SfloW2D.template" "$work_dir/"

(
    cd "$work_dir"
    python3 create_example.py 30 > generate.log
    sed -i.bak \
        -e 's/T_END=3.0D0/T_END=1.0D-3/' \
        -e 's/DT_OUTPUT=5.0D-2/DT_OUTPUT=1.0D-3/' \
        -e 's/OUTPUT_CONS_FLAG=F/OUTPUT_CONS_FLAG=T/' \
        -e 's/OUTPUT_NETCDF_FLAG=T/OUTPUT_NETCDF_FLAG=F/' \
        IMEX_SfloW2D.inp
    "$executable" > run.log
    python3 - <<'PY'
import numpy as np

state = np.loadtxt("inclinedExcavation2D_0001.q_2d")
x, y, mass = state[:, 0], state[:, 1], state[:, 2]
regions = {
    "uphill": x < 10.0,
    "south": (y < -5.0) & (x >= 10.0) & (x < 20.0),
    "north": (y >= 5.0) & (x >= 10.0) & (x < 20.0),
    "downhill": x >= 20.0,
}
outside = {name: float(np.sum(mass[mask])) for name, mask in regions.items()}
total = float(np.sum(mass))
relative_outside = sum(abs(value) for value in outside.values()) / max(total, 1.0)
if np.min(mass) < -1.0e-12 or relative_outside > 2.0e-12:
    raise SystemExit(
        f"inclined-excavation leakage: regions={outside}, "
        f"relative={relative_outside:.3e}, min_mass={np.min(mass):.3e}"
    )
print(f"outside masses [kg]: {outside}")
print(f"relative outside mass: {relative_outside:.3e}")
PY
)

echo "PASS: no resolved uphill/lateral transport across the Q1 excavation crest"
