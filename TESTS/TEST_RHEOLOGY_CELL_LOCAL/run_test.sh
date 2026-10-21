#!/bin/sh
# D02 characterization; nonzero exit means scientific acceptance is still open.
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
work_dir=${D02_WORK_DIR:-$(mktemp -d /tmp/imex-d02.XXXXXX)}
mkdir -p "$work_dir"
work_dir=$(CDPATH= cd -- "$work_dir" && pwd)
echo "D02 evidence directory: $work_dir"
compiler=${FC:-gfortran}
failed=0
for profile in strict optimized; do
    mkdir -p "$work_dir/$profile"
    cd "$work_dir/$profile"
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all -finit-real=snan' ;;
        optimized) flags='-Ofast -g -funroll-all-loops' ;;
    esac
    common='-fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -ffree-line-length-none'
    objects=''
    for source in core/parameters_2d core/diagnostics_2d numerics/lapack_interfaces_2d \
        core/constitutive_parameters_2d numerics/equation_metadata_2d core/complexify \
        geometry/geometry_2d state/state_conversion_2d physics/equation_terms_2d \
        numerics/nonlinear_solver_2d; do
        "$compiler" $flags $common -I"$repo_dir/src/state" -c "$repo_dir/src/$source.f90"
        objects="$objects ${source##*/}.o"
    done
    # IEEE checks and the test-side FD subtraction stay strict in both profiles.
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/probe_sources.f90"
    "$compiler" $common $objects probe_sources.o -llapack -o probe_sources
    if [ "$profile" = optimized ]; then
        # Distinguish an exception exposed by traps from an ordinary optimized
        # run. Production module objects are identical in both executables.
        "$compiler" -O0 -g -fcheck=all -fopenmp -ffree-line-length-none \
            -c "$test_dir/probe_sources.f90" -o probe_untrapped.o
        "$compiler" -fopenmp $objects probe_untrapped.o -llapack -o probe_untrapped
    fi
    PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/characterize.py" \
        "$work_dir/$profile/probe_sources" "$work_dir/$profile" "$profile" || failed=1
done
exit "$failed"
