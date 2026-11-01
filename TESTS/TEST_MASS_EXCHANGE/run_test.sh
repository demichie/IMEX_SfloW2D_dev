#!/bin/sh
# Physics-only mass-exchange regression. The obsolete RK/STEP executor fixture
# was retired; Lie group transfer/restart integration tests are pending.
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 -B "$test_dir/test_projection_checker.py"
python3 -B "$test_dir/check_coupled_limiters.py"
python3 -B "$test_dir/check_process_momentum.py"
echo "PASS: mass-exchange physics static checks (Lie integration coverage pending)"
